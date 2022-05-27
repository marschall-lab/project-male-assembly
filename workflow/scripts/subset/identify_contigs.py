#!/usr/bin/env python3

import pathlib as pl
import argparse as argp
import csv
import functools as fnt

import pandas as pd


Y_SPECIFIC_MOTIFS = ['DYZ1_Yq', 'DYZ18_Yq', 'DYZ3-sec_Ycentro']
UNSPECIFIC_MOTIFS = ['DYZ2_Yq']
UNSPECIFIC_THRESHOLD = 300
PRIMARY_ALN_THRESHOLD_PCT = 90


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--agg-align',
        '-a',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='align',
        required=True,
        help='Path to aggregated contig-to-reference alignment table.'
    )
    parser.add_argument(
        '--agg-motif',
        '-m',
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs='+',
        dest='motifs',
        required=True,
        help='Path to aggregated motif hits.'
    )
    parser.add_argument(
        '--out-stats',
        '-s',
        type=lambda x: pl.Path(x).resolve(),
        dest='out_stats',
        required=True,
        help='Path to chrY contig stats table output'
    )
    parser.add_argument(
        '--out-names',
        '-n',
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest='out_names',
        required=False,
        help='Path to chrY contig names (listing)'
    )
    parser.add_argument(
        '--out-bed',
        '-b',
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest='out_bed',
        required=False,
        help='Path to chrY contigs (BED format)'
    )
    parser.add_argument(
        '--select-chrom',
        '-c',
        type=str,
        default='chrY',
        choices=['chrY', 'chrX'],
        dest='chrom',
        help='Specify chromosome to select: chrY (default) or chrX'
    )
    args = parser.parse_args()
    return args


# TODO: make a generic "apply_rule" function
# and avoid code duplication
def rule_select_chrom_only_alignments(table, reason, chrom):

    chry_align_columns = [c for c in table.columns if c.startswith(f'num_aln_{chrom}')]
    other_align_columns = [c for c in table.columns if c.startswith('num_aln') and c not in chry_align_columns]

    select_y = table[chry_align_columns].sum(axis=1) > 0
    no_other = table[other_align_columns].sum(axis=1) < 1

    selector = select_y & no_other
    select, remain = apply_selector(selector, table, reason)
    return select, remain


def rule_select_y_specific_motif_hits(table, motif, reason):

    chry_align_columns = [c for c in table.columns if c.startswith('num_aln_chrY')]
    motif_hit_column = f'num_hits_hiq_{motif}'

    select_y = table[chry_align_columns].sum(axis=1) > 0
    motif_hit = table[motif_hit_column] > 0

    selector = select_y & motif_hit
    select, remain = apply_selector(selector, table, reason)
    return select, remain


def rule_select_many_unspecific_motif_hits(table, motif, reason):

    chry_primary_column = 'num_aln_chrY_PRI'
    motif_hit_column = f'num_hits_hiq_{motif}'

    select_y = table[chry_primary_column] > 0
    motif_hit = table[motif_hit_column] > UNSPECIFIC_THRESHOLD

    selector = select_y & motif_hit
    select, remain = apply_selector(selector, table, reason)
    return select, remain


def rule_select_unaligned_specific_hits(table, motif, reason):

    align_columns = [c for c in table.columns if c.startswith('num_aln_')]
    motif_hit_column = f'num_hits_hiq_{motif}'

    select_y = table[align_columns].sum(axis=1) < 1
    motif_hit = table[motif_hit_column] > 0

    selector = select_y & motif_hit
    select, remain = apply_selector(selector, table, reason)
    return select, remain


def rule_select_majority_primary_alignment(table, reason, chrom):

    chrom_primary_column = f'pct_aln_{chrom}_PRI'
    selector = table[chrom_primary_column] > PRIMARY_ALN_THRESHOLD_PCT
    select, remain = apply_selector(selector, table, reason)
    return select, remain


def apply_selector(selector, table, reason):

    if not selector.any():
        selected_contigs = None
        remaining_contigs = table
    else:
        selected_contigs = table.loc[selector, :].copy()
        remaining_contigs = table.loc[~selector, :].copy()
        selected_contigs['reason'] = reason
    return selected_contigs, remaining_contigs


def merge_aggregated_tables(align, motifs, chrom):

    keep_motif_columns = [
        'target',  # tig name
        'num_hits_hiq'
    ]

    df = pd.read_csv(align, sep='\t', header=0)
    processed_motifs = []
    if chrom == 'chrY':
        for motif_table in motifs:
            tmp = pd.read_csv(motif_table, sep='\t', header=0)
            motif_name = tmp['query'].unique()
            processed_motifs.append(motif_name)
            assert len(motif_name) == 1
            motif_name = motif_name[0]
            drop_columns = [c for c in tmp.columns if c not in keep_motif_columns]
            tmp.drop(drop_columns, axis=1, inplace=True)
            tmp.rename(
                {
                    'target': 'tig',
                    'num_hits_hiq': f'num_hits_hiq_{motif_name}'
                },
                axis=1,
                inplace=True
            )
            df = df.merge(tmp, how='outer', left_on='tig', right_on='tig')
        df.fillna(0, inplace=True)
    return df, processed_motifs


def select_chrom_contigs(table, chrom):

    selected = []
    remaining = table.copy()

    reason = f"Rule 1: only {chrom} alignments"
    subset, remaining = rule_select_chrom_only_alignments(
        remaining,
        reason,
        chrom
    )
    if subset is not None:
        print(reason, ': ', subset.shape[0])
        selected.append(subset)

    if chrom == 'chrY':

        for motif in Y_SPECIFIC_MOTIFS:
            reason = f"Rule 2: mixed alignments with specific motif hits [{motif}]"
            subset, remaining = rule_select_y_specific_motif_hits(
                remaining,
                motif,
                reason
            )
            if subset is not None:
                print(reason, ': ', subset.shape[0])
                selected.append(subset)

        for motif in UNSPECIFIC_MOTIFS:
            reason = f"Rule 3: mixed alignments (primary to chrY) with many (>{UNSPECIFIC_THRESHOLD}) unspecific motif hits [{motif}]"
            subset, remaining = rule_select_many_unspecific_motif_hits(
                remaining,
                motif,
                reason
            )
            if subset is not None:
                print(reason, ': ', subset.shape[0])
                selected.append(subset)

        for motif in Y_SPECIFIC_MOTIFS:
            reason = f"Rule 4: unaligned contigs with specific motif hits [{motif}]"
            subset, remaining = rule_select_unaligned_specific_hits(
                remaining,
                motif,
                reason
            )
            if subset is not None:
                print(reason, ': ', subset.shape[0])
                selected.append(subset)

    reason = f"Rule 5: majority (>{PRIMARY_ALN_THRESHOLD_PCT}% bp) of contig primary alignments to {chrom}"
    subset, remaining = rule_select_majority_primary_alignment(remaining, reason, chrom)
    if subset is not None:
        print(reason, ': ', subset.shape[0])
        selected.append(subset)

    if not selected:
        raise RuntimeError(f'No {chrom} contigs selected')
    selected = pd.concat(selected, axis=0, ignore_index=False)
    selected['reason'] = '"' + selected['reason'] + '"'
    num_columns = [c for c in selected.columns if c not in ['tig', 'reason']]
    selected[num_columns] = selected[num_columns].astype(int)
    return selected


def keep_stats_columns(chrom, column):

    if column in ['tig', 'tig_length', 'reason']:
        keep = True
    elif 'hits' in column:
        keep = True
    elif chrom in column:
        keep = True
    else:
        keep = False
    return keep


def main():

    args = parse_command_line()
    df, motifs = merge_aggregated_tables(args.align, args.motifs, args.chrom)
    selected = select_chrom_contigs(df, args.chrom)

    keep_chrom_stats_columns = fnt.partial(keep_stats_columns, args.chrom)
    keep_columns = list(filter(keep_chrom_stats_columns, selected.columns))

    args.out_stats.parent.mkdir(parents=True, exist_ok=True)
    selected.sort_values(['tig_length', 'tig'], inplace=True, ascending=False)
    selected.to_csv(
        args.out_stats,
        sep='\t',
        header=True,
        index=False,
        columns=keep_columns,
        quoting=csv.QUOTE_NONE
    )

    if args.out_names is not None:
        args.out_names.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_names, 'w') as dump:
            _ = dump.write('\n'.join(sorted(selected['tig'])))

    if args.out_bed is not None:
        args.out_bed.parent.mkdir(parents=True, exist_ok=True)
        bed_regions = sorted([f"{tig}\t0\t{int(tig_length)}" for tig, tig_length in selected[['tig', 'tig_length']].itertuples(index=False)])
        with open(args.out_bed, 'w') as dump:
            _ = dump.write('\n'.join(bed_regions) + '\n')

    return 0


if __name__ == '__main__':
    main()
