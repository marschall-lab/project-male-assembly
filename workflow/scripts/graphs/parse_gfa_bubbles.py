#!/usr/bin/env python3

import argparse as argp
import csv as csv
import hashlib as hl
import pathlib as pl
import io
import sys

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--bubble-bed",
        "--bubbles",
        "-bed",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="bubbles",
        help="GFAtools bubbles output"
    )
    parser.add_argument(
        "--node-annotation",
        "-na",
        "-a",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="annotation",
        help="CSV file with node annotations (for Bandage)"
    )
    parser.add_argument(
        "--out-table",
        "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_table",
        help="Path to output file (TSV)"
    )
    parser.add_argument(
        "--out-path-seq",
        "-seq",
        "-s",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_seqs",
        help="Path to output file for path seqs (FASTA)"
    )
    args = parser.parse_args()
    return args


def parse_bubble_table(bubbles_file):

    # to fit long node sequences
    csv.field_size_limit(sys.maxsize)

    table_header = [
        "chrom", "start", "end",
        "all_nodes", "all_paths", "has_inversion",
        "shortest_path_bp", "longest_path_bp",
        "drop1", "drop2", "drop3",
        "nodelist", "shortest_path_seq", "longest_path_seq"
    ]

    seq_buffer = io.StringIO()

    bubbles = []
    with open(bubbles_file, "r", newline="") as table:
        reader = csv.DictReader(table, fieldnames=table_header, delimiter="\t")
        for ln, row in enumerate(reader, start=1):
            bubble_num = ln
            bubble_alleles = int(row["all_nodes"]) - 2
            assert bubble_alleles > 0, row
            short_length = 0 if row["shortest_path_bp"] is None else int(row["shortest_path_bp"])
            long_length = 0 if row["longest_path_bp"] is None else int(row["longest_path_bp"])
            bubble_midsize = (short_length + long_length) // 2
            try:
                bubble_id = hl.md5(row["nodelist"].encode("utf-8")).hexdigest()
            except AttributeError:
                print(row)
                raise
            if row["shortest_path_seq"] != "*":
                seq_buffer.write(f">{bubble_num}_shortest_{bubble_id}\n")
                seq_buffer.write(row["shortest_path_seq"] + "\n")
            seq_buffer.write(f">{bubble_num}_longest_{bubble_id}\n")
            seq_buffer.write(row["longest_path_seq"] + "\n")
            for pos, node in enumerate(row["nodelist"].split(","), start=0):
                node_info = {
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "node": node,
                    "is_source": 0,
                    "is_sink": 0,
                    "bubble_num": bubble_num,
                    "bubble_id": bubble_id,
                    "bubble_midsize": bubble_midsize,
                    "bubble_has_inversion": int(row["has_inversion"]),
                    "bubble_alleles": bubble_alleles
                }
                if pos == 0:
                    node_info["is_source"] = 1
                elif pos == int(row["all_nodes"]) - 1:
                    node_info["is_sink"] = 1
                else:
                    pass
                bubbles.append(node_info)
    df = pd.DataFrame.from_records(bubbles)
    return df, seq_buffer


def main():
    args = parse_command_line()
    annotation = pd.read_csv(args.annotation, header=0)
    # this is only properly defined for T2T-Y
    annotation.drop(["ref_start", "ref_end"], axis=1, inplace=True)
    bubbles, path_seqs = parse_bubble_table(args.bubbles)
    bubbles = bubbles.merge(annotation, on="node", how="outer")
    
    args.out_table.parent.mkdir(exist_ok=True, parents=True)
    bubbles.to_csv(args.out_table, sep="\t", header=True, index=False)

    args.out_seqs.parent.mkdir(exist_ok=True, parents=True)
    with open(args.out_seqs, "w") as fasta:
        _ = fasta.write(path_seqs.getvalue())
    
    return 0


if __name__ == "__main__":
    main()