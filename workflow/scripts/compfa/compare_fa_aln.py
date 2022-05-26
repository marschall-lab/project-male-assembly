#!/usr/bin/env python3

import os
import sys
import pathlib as pl
import argparse as argp
import itertools as itt
import io
import gzip


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--iput',
        '-i',
        dest='input',
        nargs='+',
        type=lambda x: pl.Path(x).resolve(strict=True),
        help='Path to input FASTA files (must have extension .gz)'
    )
    parser.add_argument(
        '--show-id-pos',
        '-s',
        action='store_true',
        default=False,
        dest='show_id_pos',
        help='Include positions in output where all sequences are identical'
    )
    parser.add_argument(
        '--output',
        '-o',
        nargs='?',
        type=lambda x: pl.Path(x).resolve(),
        default=sys.stdout,
        help='Specify output. Default: stdout'
    )
    parser.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        dest='quiet',
        help='Do not print debug messages to stderr'
    )
    parser.add_argument(
        '--first-index',
        '-f',
        type=int,
        default=1,
        dest='first_index',
        help='Index of first position in enumeration, e.g., 0 or 1. Default: 1'
    )
    parser.add_argument(
        '--field-sep',
        '-d',
        type=str,
        default='\t',
        dest='field_sep',
        help='Specify output field separator. Default: <TAB>'
    )
    parser.add_argument(
        '--uppercase',
        '-u',
        dest='uppercase',
        action='store_true',
        default=False,
        help='Convert FASTA sequences to uppercase symbols. Default: False'
    )

    args = parser.parse_args()
    return args


def determine_opener(filepath, quiet):
    file_ext = filepath.suffix
    if file_ext in ['.gz', '.gzip']:
        msg = f'Detected gzip-compressed input: {filepath}'
        if not quiet:
            sys.stderr.write(f'\n{msg}\n')
        open_fun, open_mode = gzip.open, 'rt'
    elif file_ext in ['.fa', '.fasta', 'fna']:
        open_fun, open_mode = open, 'r'
    else:
        msg = f'Non-standard FASTA file extension: {filepath}'
        if not quiet:
            sys.stderr.write(f'\n{msg}\n')
        open_fun, open_mode = open, 'r'
    return open_fun, open_mode


def load_fasta_sequences(input_files, uppercase, quiet):

    fasta_seqs = []
    seq_names = []
    seq_lengths = []
    for filepath in input_files:
        open_fun, open_mode = determine_opener(filepath, quiet)
        current_seq = ''
        current_name = None
        with open_fun(filepath, open_mode) as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if current_name is not None:
                        seq_names.append(current_name)
                        seq_lengths.append(len(current_seq))
                        if uppercase:
                            fasta_seqs.append(current_seq.upper())
                        else:
                            fasta_seqs.append(current_seq)
                    current_name = line.strip().split()[0][1:]
                    current_seq = ''
                else:
                    current_seq += line.strip()
            if current_seq:
                seq_names.append(current_name)
                seq_lengths.append(len(current_seq))
                if uppercase:
                    fasta_seqs.append(current_seq.upper())
                else:
                    fasta_seqs.append(current_seq)

    if len(set(seq_lengths)) > 1:
        msg = 'Warning: FASTA sequences have different lengths'
        sys.stderr.write(f'\n{msg}\n')
        for sn, sl in zip(seq_names, seq_lengths):
            sys.stderr.write(f'{sn}\t{sl}\n')

    assert len(seq_names) == len(seq_lengths) == len(fasta_seqs)

    return fasta_seqs, seq_names


def enum_all(sequences, first_index):

    for pos, values in enumerate(itt.zip_longest(*tuple(sequences), fillvalue='-'), start=first_index):
        yield pos, values
    return


def enum_different(sequences, first_index):

    for pos, values in enumerate(itt.zip_longest(*tuple(sequences), fillvalue='-'), start=first_index):
        if len(set(values)) == 1:
            continue
        yield pos, values
    return


def main():
    args = parse_command_line()

    fasta_seqs, seq_names = load_fasta_sequences(args.input, args.uppercase, args.quiet)

    iter = enum_all if args.show_id_pos else enum_different

    fs = args.field_sep

    if not isinstance(args.output, io.TextIOWrapper):
        args.output.parent.mkdir(parents=True, exist_ok=True)

    header = f'#POS{fs}' + f'{fs}'.join(seq_names) + '\n'
    try:
        with open(args.output, 'w') as dump:
            dump.write(header)
            for pos, values in iter(fasta_seqs, args.first_index):
                dump.write(f'{pos}{fs}' + f'{fs}'.join(values) + '\n')
    except TypeError:
        with args.output as dump:
            for pos, values in iter(fasta_seqs, args.first_index):
                dump.write(f'{pos}{fs}' + f'{fs}'.join(values) + '\n')

    return


if __name__ == '__main__':
    main()
