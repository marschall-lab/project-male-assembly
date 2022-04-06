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
        required=True,
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
        required=True,
        help='Path to FASTA output file to dump adapted subset (chromosome) assembly.'
    )
    args = parser.parse_args()
    return args


def main():

    args = parse_command_line()

    name_map = json.load(open(args.names, 'r'))

    revcomp_table = str.maketrans(COMPLEMENT_NUCLEOTIDES)

    args.out_wg.parent.mkdir(parents=True, exist_ok=True)

    subset_buffer = []
    with dnaio.FastaWriter(args.out_wg) as fasta_wg:
        with dnaio.open(args.input_fasta) as fasta_in:
            for record in fasta_in:
                if NUCLEOTIDES.match(record.sequence) is None:
                    nucs = col.Counter(record.sequence.upper())
                    raise ValueError(f'Illegel characters in sequence record {record.name}: {nucs}')
                out_name = name_map.get(record.name, record.name)
                if '.RV.' in out_name:
                    out_seq = record.sequence.translate(revcomp_table)[::-1]
                else:
                    out_seq = record.sequence
                fasta_wg.write(out_name, out_seq)
                if out_name != record.name:
                    subset_buffer.append((out_name, out_seq))

    args.out_sub.parent.mkdir(parents=True, exist_ok=True)
    with dnaio.FastaWriter(args.out_sub) as fasta_sub:
        for name, seq in sorted(subset_buffer):
            fasta_sub.write(name, seq)

    return 0


if __name__ == '__main__':
    main()
