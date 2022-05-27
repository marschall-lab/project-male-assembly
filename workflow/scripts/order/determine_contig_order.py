#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import collections as col
import re
import json

import pandas as pd


SAMPLE_NAME = re.compile('[A-Z0-9a-z]+')


PAF_INPUT_COLUMNS = [
    'query',
    'query_length',
    'query_aln_start',
    'query_aln_end',
    'orientation',
    'target',
    'target_length',
    'target_aln_start',
    'target_aln_end',
    'sequence_matches',
    'aln_block_length',
    'mapq',
    'tag_NM',
    'tag_ms',
    'tag_AS',
    'tag_nn',
    'tag_tp',
    'tag_cm',
    'tag_s1',
    'tag_s2',
    'tag_de',
    'tag_rl',
    'tag_cs',
    'tag_cg'
]

KEEP_PAF_COLUMNS = [
    'query',
    'query_length',
    'query_aln_start',
    'query_aln_end',
    'target',
    'target_length',
    'target_aln_start',
    'target_aln_end',
    'sequence_matches',
    'tag_tp'
]
assert all(c in PAF_INPUT_COLUMNS for c in KEEP_PAF_COLUMNS)


COVERAGE_INPUT_COLUMNS = [
    'region',
    'region_start',
    'region_end',
    'region_name',
    'target',
    'target_start',
    'target_end',
    'query',
    'mapq',
    'orientation',
    'overlap'
]

KEEP_COVERAGE_COLUMNS = [
    'region_name',
    'target',
    'target_start',
    'target_end',
    'query',
    'mapq',
    'orientation',
    'overlap'
]
assert all(c in COVERAGE_INPUT_COLUMNS for c in KEEP_COVERAGE_COLUMNS)


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--names',
        '-n',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='names',
        required=True,
        help='Path to text file listing identified chrY contig names.'
    )
    parser.add_argument(
        '--class-coverage',
        '-c',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='coverage',
        required=True,
        help='Path to TSV file containing intersect beween Y seq. classes and contig alignments.'
    )
    parser.add_argument(
        '--seq-classes',
        '-s',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='seq_classes',
        required=True,
        help='Path to annotation file containing Y sequence classes.'
    )
    parser.add_argument(
        '--paf',
        '-p',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='paf',
        required=True,
        help='Path to PAF contig-to-reference alignment file (can be gzipped).'
    )
    parser.add_argument(
        '--sample-name',
        '-sn',
        type=str,
        dest='sample_name',
        required=True,
        help='(Well-bahaved) sample name to be appended to the new contig names'
    )
    parser.add_argument(
        '--min-cov-threshold',
        '-t',
        type=int,
        default=95,
        dest='min_cov_threshold',
        help='Specify minimal contig coverage threshold to select regions, '
             'i.e. chrY seq. classes, for contig renaming. Default: 95 (%%)'
    )
    parser.add_argument(
        '--output',
        '-o',
        type=lambda x: pl.Path(x),
        required=True,
        help='Path to output file listing new contig names.'
    )
    parser.add_argument(
        '--dump-mappings',
        '-map',
        type=str,
        nargs='*',
        default=None,
        dest='dump_maps',
        choices=['tsv', 'TSV', 'json', 'JSON', 'sed', 'SED'],
        help='If specified, dump files mapping old-to-new/new-to-old contig names. '
             'This will strip off the file extension of "--output" and use that as '
             'base name with extension "otn-map.[EXT]" and "nto-map.[EXT]".'
    )
    parser.add_argument(
        '--process-chrom',
        '-chrom',
        type=str,
        default='chrY',
        choices=['chrY', 'chrX'],
        dest='chrom',
        help='Specify which chromosome to process: chrY (default) or chrX'
    )
    args = parser.parse_args()
    return args


def load_sequence_class_order(file_path):
    """
    Loads Y seq. class annoation file
    to set global order for regions
    """
    name_to_ordernum = dict()
    with open(file_path, 'r') as bed_file:
        for order, line in enumerate(bed_file, start=1):
            order_num = f'{order:02}'
            name = line.split()[3]
            name_to_ordernum[name] = order_num
    return name_to_ordernum


def load_sequence_class_coverage(file_path, y_contigs, seq_class_to_order):
    """
    This loads the output TSV of rule

    intersect_contig_align_chry_seq_classes

    which intersects (primary) contig alignments
    with the seq. class annotation file, but on the
    whole-genome level; hence,
    1. select only target chrY
    2. select only query "y contigs" (as determined before)
    """

    df = pd.read_csv(
        file_path, sep='\t', header=None, index_col=False,
        names=COVERAGE_INPUT_COLUMNS, usecols=KEEP_COVERAGE_COLUMNS
    )
    df = df.loc[(df['target'] == 'chrY') & (df['query'].isin(y_contigs)), :].copy()
    df['orientation'].replace({'+': 'FW', '-': 'RV'}, inplace=True)
    df['order'] = df['region_name'].replace(seq_class_to_order)
    df['order_num'] = df['order'].astype(int)
    return df


def load_alignment_file(file_path, y_contigs):
    """
    """
    df = pd.read_csv(
        file_path, sep='\t',
        index_col=False, header=None,
        names=PAF_INPUT_COLUMNS, usecols=KEEP_PAF_COLUMNS
    )
    # drop everything not aligned to chrY
    df.drop(df.index[df['target'] != 'chrY'], axis=0, inplace=True)
    # consider only primary alignments
    df.drop(df.index[~df['tag_tp'].isin(['tp:A:P', 'tp:A:p'])], axis=0, inplace=True)
    df.drop('tag_tp', inplace=True, axis=1)

    paf_contigs = set(df['query'].unique())
    aligned_contigs = paf_contigs.intersection(y_contigs)
    unaligned_contigs = set(y_contigs) - paf_contigs

    return df, aligned_contigs, unaligned_contigs


def determine_relevant_regions(tig, tig_length, region_coverages, min_cov_t):
    """
    2022-05-27
    Exclude regions hit by spurious alignments that are non-neighboring
    relative to the bulk of the aligned tig sequence
    """
    region_coverages['pct_cov_query'] = (region_coverages['overlap'] / tig_length * 100).round(1)
    region_pct_cov = region_coverages.groupby('order_num')['pct_cov_query'].sum()

    all_order_nums = sorted(set(region_pct_cov.index))

    # identify all subsequences of consecutive regions (chrY seq. classes)
    subseq_num = 0
    subsequences = col.defaultdict(list)
    for on in all_order_nums:
        if len(all_order_nums) == 1:
            subseq_num += 1
            subsequences[subseq_num].append(on)
        elif on - 1 in all_order_nums:
            # not a starting region
            subsequences[subseq_num].append(on)
        elif on - 1 in all_order_nums and on + 1 in all_order_nums:
            # region in the middle is not a starting region
            subsequences[subseq_num].append(on)
        elif on + 1 in all_order_nums:
            # starting region
            subseq_num += 1
            subsequences[subseq_num].append(on)
        else:
            # isolated region, no predecessor or successor,
            # is a starting region
            subseq_num += 1
            subsequences[subseq_num].append(on)

    # sum up pct. coverage per subsequences
    # threshold set at >95% to select subsequence(s)
    # used for contig renaming
    cov_per_subsequence = []
    for _, order_nums in subsequences.items():
        total = region_pct_cov.loc[order_nums].sum()
        cov_per_subsequence.append((total, order_nums))

    # now select all subsequences needed to get over 95% contig coverage
    selected_order_nums = []
    cumsum = 0
    for cov, order_nums in sorted(cov_per_subsequence, reverse=True):
        cumsum += cov
        selected_order_nums.extend(order_nums)
        if cumsum > min_cov_t:
            break

    return selected_order_nums


def find_contig_order(seq_class_coverage, alignments, sample_name, min_cov_t):

    renamed = col.defaultdict(set)
    for tig, coverages in seq_class_coverage.groupby('query'):
        # update 2022-05-27
        # compute pct overlap for query (assembled tig) to determine
        # if certain regions (alignments) should not be considered for
        # renaming because they are distant relative to the others, and
        # cover only a small percentage of the query contig
        query_length = alignments.loc[alignments['query'] == tig, 'query_length'].values[0]
        selected_order_nums = determine_relevant_regions(tig, query_length, coverages.copy(), min_cov_t)

        # reduce considered coverages table to only selected regions
        relevant_coverages = coverages.loc[coverages['order_num'].isin(selected_order_nums), :]

        # determine T2T-relative orientation as simple majority
        # update 2022-05-27: this may result in undesired revcomp if the assembled
        # tig happens to coincide with an actual inversion; these case (if report
        # by Pille H.) are hard-coded to prevent a revcomp operation
        majority_orientation = relevant_coverages.groupby('orientation')['overlap'].sum().idxmax()

        # lowest start in alignment file: needed to order tigs
        # that span the same region(s), which would otherwise
        # not lead to a proper sort order
        # lowest_start = alignments.loc[alignments['query'] == tig, 'target_aln_start'].min()
        # update 2022-05-27
        # since the PAF alignment file does not contain the info about the
        # chrY sequence classes, replace "lowest start" (above) by "lowest start"
        # in relevant_coverages. That could be more likely to create identical
        # aln. starts
        lowest_start = relevant_coverages['target_start'].min()

        uniq = coverages.drop_duplicates(['region_name', 'order_num'])

        uniq = relevant_coverages.drop_duplicates(['region_name', 'order_num'])

        # order_num = int(order), to make min/max work on numbers
        lowest = uniq['order_num'].idxmin()
        highest = uniq['order_num'].idxmax()
        cls_start = uniq.at[lowest, 'region_name']
        cls_end = uniq.at[highest, 'region_name']
        order_start = uniq.at[lowest, 'order']
        order_end = uniq.at[highest, 'order']
        order_prefix = f'chrY.{order_start}-{order_end}'
        class_suffix = f'{cls_start}-{cls_end}.{majority_orientation}.{tig}.{sample_name}'

        # lowest start: smallest alignment coordinate, this is used
        # to sort all tigs spanning the same region(s)
        renamed[order_prefix].add((lowest_start, class_suffix, tig))

    # now derive final new names
    new_names = []
    for prefix, suffices in renamed.items():
        for pos, (_, suffix, tig) in enumerate(sorted(suffices), start=1):
            new_name = f'{prefix}.{pos:02}.{suffix}'
            new_names.append((new_name, tig))

    assert new_names
    return new_names


def enumerate_chrx_contigs(paf_file, contigs, sample_name):

    # a bit of code duplication here...
    df = pd.read_csv(
        paf_file, sep='\t',
        index_col=False, header=None,
        names=PAF_INPUT_COLUMNS, usecols=KEEP_PAF_COLUMNS + ['orientation']
    )
    # drop everything not aligned to chrY
    df.drop(df.index[df['target'] != 'chrX'], axis=0, inplace=True)
    # consider only primary alignments
    df.drop(df.index[~df['tag_tp'].isin(['tp:A:P', 'tp:A:p'])], axis=0, inplace=True)
    df.drop('tag_tp', inplace=True, axis=1)
    df['target_aln_length'] = (df['target_aln_end'] - df['target_aln_start']).astype(int)
    df['orientation'].replace({'+': 'FW', '-': 'RV'}, inplace=True)

    name_parts = []
    for tig_name, alignments in df.groupby('query'):
        if tig_name not in contigs:
            # this is important: contigs can only be identified
            # as chrX based on the alignment, but not all contigs
            # having some (small...) alignment to chrX are also
            # chrX contigs (that decision is made in "identify_contigs.py").
            # Hence, process only what is in the input here.
            continue
        majority_orientation = alignments.groupby('orientation')['target_aln_length'].sum().idxmax()
        lowest_start = alignments['target_aln_start'].min()
        name_parts.append((lowest_start, majority_orientation, tig_name))

    new_names = []
    renamed = set()
    for sort_order, (ls, orient, tig_name) in enumerate(sorted(name_parts), start=1):
        if len(name_parts) < 100:
            new_name = f'chrX.{sort_order:02}.{orient}.{tig_name}.{sample_name}'
        else:
            new_name = f'chrX.{sort_order:03}.{orient}.{tig_name}.{sample_name}'
        new_names.append((new_name, tig_name))
        renamed.add(tig_name)

    assert new_names
    return new_names


def dump_output_files(new_contig_names, output_file, output_maps):

    output_maps = [m.lower() for m in output_maps]

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as dump:
        _ = dump.write('\n'.join([t[0] for t in new_contig_names]) + '\n')

    if 'tsv' in output_maps:
        otn_file = output_file.with_suffix('.otn-map.tsv')
        nto_file = output_file.with_suffix('.nto-map.tsv')
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f'{old_name}\t{new_name}\n')
                    _ = nto.write(f'{new_name}\t{old_name}\n')

    if 'json' in output_maps:
        otn_file = output_file.with_suffix('.otn-map.json')
        nto_file = output_file.with_suffix('.nto-map.json')

        otn_map = dict((old_name, new_name) for new_name, old_name in new_contig_names)
        with open(otn_file, 'w') as otn:
            json.dump(otn_map, otn, indent=1)

        nto_map = dict(new_contig_names)
        with open(nto_file, 'w') as nto:
            json.dump(nto_map, nto, indent=1)

    if 'sed' in output_maps:
        otn_file = output_file.with_suffix('.otn-map.sed')
        nto_file = output_file.with_suffix('.nto-map.sed')
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f's/\\b{old_name}\\b/{new_name}/g\n')
                    _ = nto.write(f's/\\b{new_name}\\b/{old_name}/g\n')

    return


def main():

    args = parse_command_line()

    # check that sample name is not sth crazy
    if SAMPLE_NAME.match(args.sample_name) is None:
        raise ValueError(f'Sample name should only be letters and numbers: {args.sample_name}')

    # original contig names
    chrom_contigs = open(args.names).read().strip().split()
    new_contig_names = []

    if args.chrom == 'chrY':
        # load global order
        seq_class_to_order = load_sequence_class_order(args.seq_classes)

        # contig coverage in Y seq. classes / regions
        seq_class_coverage = load_sequence_class_coverage(
            args.coverage,
            chrom_contigs,
            seq_class_to_order
        )

        alignments, chrom_contigs, unaligned_contigs = load_alignment_file(args.paf, chrom_contigs)

        if unaligned_contigs:
            prefix = 'chrY.99-99.00.UNK-UNK.NA.'
            for contig in unaligned_contigs:
                new_name = prefix + contig + f'.{args.sample_name}'
                new_contig_names.append((new_name, contig))

        new_names_aln = find_contig_order(seq_class_coverage, alignments, args.sample_name, args.min_cov_threshold)

    elif args.chrom == 'chrX':
        # chrX is much simpler because no sequence class annotation is available,
        # and, by construction, no unaligned contigs can be identified
        new_names_aln = enumerate_chrx_contigs(args.paf, chrom_contigs, args.sample_name)

    else:
        raise ValueError(f'Unknown chromosome: {args.chrom}')

    new_contig_names = sorted(new_contig_names + new_names_aln)

    dump_output_files(
        new_contig_names,
        args.output,
        [] if args.dump_maps is None else args.dump_maps
    )

    return 0


if __name__ == '__main__':
    main()
