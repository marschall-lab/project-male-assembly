#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import collections as col

import pandas as pd


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

KEEP_COLUMNS = [
    'query',
    'query_length',
    'query_aln_start',
    'query_aln_end',
    'target',
    'target_length',
    'target_aln_start',
    'target_aln_end',
    'sequence_matches',
    'mapq',
    'tag_AS',
    'tag_tp'
]
assert all(c in PAF_INPUT_COLUMNS for c in KEEP_COLUMNS)


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--paf',
        '-p',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='paf',
        help='Path to PAF contig-to-reference alignment file (can be gzipped).'
    )
    parser.add_argument(
        '--output',
        '-o',
        type=lambda x: pl.Path(x),
        help='Path to output table (tab-separated).'
    )
    args = parser.parse_args()
    return args


def norm_as_tag(value):
    """
    Value may be nan for unaligned contigs
    """
    try:
        i = int(value.split(':')[-1])
    except AttributeError:
        i = 0
    return i


def get_tp_tag_mapping():
    lut = {
        'tp:A:P': 'PRI',
        'tp:A:p': 'PRI',
        'tp:A:S': 'SEC',
        'tp:A:s': 'SEC',
        'tp:A:I': 'INV',
        'tp:A:i': 'INV'
    }
    return lut


def group_chromosomes(value):
    if value == 'chrY':
        return 'chrY'
    elif value == 'chrX':
        return 'chrX'
    elif value == '*':
        return 'UNK'
    else:
        return 'chrN'


def main():

    args = parse_command_line()

    ctg_aln = pd.read_csv(args.paf, sep='\t', header=None, names=PAF_INPUT_COLUMNS, usecols=KEEP_COLUMNS, index_col=False)
    ctg_aln['aln_score'] = ctg_aln['tag_AS'].apply(norm_as_tag)
    tp_mapping = get_tp_tag_mapping()
    ctg_aln['aln_type'] = ctg_aln['tag_tp'].apply(lambda x: tp_mapping.get(x, 'UNK'))
    ctg_aln['chrom_group'] = ctg_aln['target'].apply(group_chromosomes)

    ctg_aln.drop(['tag_AS', 'tag_tp', 'target'], axis=1, inplace=True)

    # Aggregate alignment records per grouping:
    # query = assembled utig name
    # group = reference chromosome group
    # aln_type = alignment type (primary, secondary etc.)
    agg_records = col.defaultdict(dict)
    for (query, group, aln_type), alignments in ctg_aln.groupby(['query', 'chrom_group', 'aln_type']):
        query_length = alignments.at[alignments.index[0], 'query_length']
        if group == 'UNK':
            # unaligned contig, all stats zero
            agg_records[query].update({'tig': query, 'tig_length': query_length})
            continue
        aln_bp = (alignments['query_aln_end'] - alignments['query_aln_start']).sum()
        aln_pct = round(aln_bp / query_length * 100, 0)
        record = {
            'tig': query,
            'tig_length': query_length,
            f'num_aln_{group}_{aln_type}': alignments.shape[0],
            f'pct_aln_{group}_{aln_type}': aln_pct,
            f'mean_mapq_{group}_{aln_type}': alignments['mapq'].mean(),
            f'median_mapq_{group}_{aln_type}': alignments['mapq'].median(),
        }
        agg_records[query].update(record)

    agg_records = pd.DataFrame.from_records(list(agg_records.values()))
    # set all missing stats to zero
    agg_records.fillna(0, inplace=True)
    agg_records.sort_values('tig', inplace=True)
    num_columns = [c for c in agg_records.columns if c != 'tig']
    agg_records[num_columns] = agg_records[num_columns].astype(int)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    agg_records.to_csv(args.output, sep='\t', header=True, index=False)
    return 0


if __name__ == '__main__':
    main()
