#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import collections as col
import re
import json

import dnaio


NUCLEOTIDES = re.compile('^[ACGT]+$', flags=re.IGNORECASE)

COMPLEMENT_NUCLEOTIDES = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'a': 't',
    'c': 'g',
    'g': 'c',
    't': 'a'
}


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--name-map',
        '-n',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='names',
        nargs='+',
        help='Path to JSON file containing chrY contig name mapping (old to new).'
    )
    parser.add_argument(
        '--input-fasta',
        '-i',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='input_fasta',
        required=True,
        help='Path to FASTA file containing unprocessed Verkko whole-genome assembly.'
    )
    parser.add_argument(
        '--out-wg',
        '-w',
        type=lambda x: pl.Path(x).resolve(),
        dest='out_wg',
        required=True,
        help='Path to FASTA output file to dump adapted whole-genome assembly.'
    )
    parser.add_argument(
        '--out-sub',
        '-s',
        type=lambda x: pl.Path(x).resolve(),
        dest='out_sub',
        nargs='+',
        help='Path to FASTA output file to dump adapted subset (chromosome) assembly.'
    )
    args = parser.parse_args()
    return args


def main():

    args = parse_command_line()
    if len(args.names) != len(args.out_sub):
        raise ValueError(f'Must specify one name mapping JSON per subset assembly FASTA:\n{args.names}\n{args.out_sub}')

    name_maps = dict()
    subset_names = dict()
    for map_idx, name_map in enumerate(args.names):
        name_map = json.load(open(name_map, 'r'))
        collisions = set(name_maps.keys()).intersection(set(name_map.keys()))
        if len(collisions) > 0:
            raise RuntimeError(f'Name collision for chromosome subsets: {collisions} / {args.input_fasta}')
        name_maps.update(name_map)
        subset_names.update(dict((new_name, map_idx) for new_name in name_map.values()))

    revcomp_table = str.maketrans(COMPLEMENT_NUCLEOTIDES)

    args.out_wg.parent.mkdir(parents=True, exist_ok=True)

    subset_buffers = col.defaultdict(list)
    with dnaio.FastaWriter(args.out_wg) as fasta_wg:
        with dnaio.open(args.input_fasta) as fasta_in:
            for record in fasta_in:
                if NUCLEOTIDES.match(record.sequence) is None:
                    nucs = col.Counter(record.sequence.upper())
                    raise ValueError(f'Illegel characters in sequence record {record.name}: {nucs}')
                out_name = name_maps.get(record.name, record.name)

                # hard-coded exception (see gh#3)
                # this contig represents an inversion relative to T2T-Y,
                # and should thus be reversed despite aligning in forward
                # orientation (request by Pille H.)
                if 'unassigned-0009118.HG03492' in out_name:
                    if '.FW.':
                        out_seq = record.sequence.translate(revcomp_table)[::-1]
                        out_name = out_name.replace('.FW.', '.RV.')
                    else:
                        # hm, quite unexpected, but should be as intended then
                        pass
                elif '.RV.' in out_name:
                    out_seq = record.sequence.translate(revcomp_table)[::-1]               
                else:
                    out_seq = record.sequence
                fasta_wg.write(out_name, out_seq)
                if out_name != record.name:
                    subset_idx = subset_names[out_name]
                    subset_buffers[subset_idx].append((out_name, out_seq))

    for subset_idx, output_subset in enumerate(args.out_sub):
        output_subset.parent.mkdir(parents=True, exist_ok=True)
        with dnaio.FastaWriter(output_subset) as fasta_sub:
            print(output_subset)
            for name, seq in sorted(subset_buffers[subset_idx]):
                fasta_sub.write(name, seq)

    return 0


if __name__ == '__main__':
    main()
