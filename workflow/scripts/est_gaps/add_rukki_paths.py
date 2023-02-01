#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import hashlib as hl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--unitig-table",
        "-u",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="unitig_table",
        help="Path to unitig table (output of 'merge_graph_infos.py')"
    )

    parser.add_argument(
        "--rukki-paths",
        "-p",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="rukki_paths",
        help="Path to Rukki paths output file."
    )

    parser.add_argument(
        "--rukki-default-gap",
        "-d",
        type=int,
        default=1000,
        dest="default_gap",
        help="Set default gap length used in Rukki. Default: 1000"
    )

    parser.add_argument(
        "--out-table",
        "-o",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_table",
        help="Path to output file for merged table."
    )

    parser.add_argument(
        "--summary-subset",
        "-sub",
        type=str,
        nargs="+",
        default=["chrY"],
        dest="summary_subset",
        help="Specify subset for dumping summary info. Default: chrY"
    )

    parser.add_argument(
        "--out-summary",
        "-s",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_summary",
        help="Path to output summary file."
    )
    args = parser.parse_args()
    return args


def parse_rukki_path(default_gap_length, unitig_infos, path_row):

    path_info = {
        "path_id": None,
        "assigned": None,
        "path_type": "disconnected",
        "contig_path": None,
        "nodes_num": 0,
        "total_gaps_num": 0,
        "total_gap_length_hpc": 0,
        "default_gaps_num": 0,
        "default_gap_length_hpc": 0,
        "est_gaps_num": 0,
        "est_gap_length_hpc": 0,
        "total_path_length_hpc": 0,
        "gapless_path_length_hpc": 0,
        "gapless_path_length_pct": 0,
    }

    if path_row.assignment == "MATERNAL":
        path_info["assigned"] = "chrX"
    elif path_row.assignment == "PATERNAL":
        path_info["assigned"] = "chrY"
    else:
        path_info["assigned"] = "UNK"
    
    path_id = hl.md5(path_row.path.encode("utf-8")).hexdigest()
    path_id = f"PATH_{path_id[:8]}"

    path_info["path_id"] = path_id

    path_members = set()
    contig_path = []
    for component in path_row.path.split(","):
        member = component.strip("+-")
        if member.startswith("["):
            gap_length = int(member.strip("[]N"))
            contig_path.append("GAP")
            path_info["total_gaps_num"] += 1
            path_info["total_gap_length_hpc"] += gap_length
            path_info["total_path_length_hpc"] += gap_length
            if gap_length == default_gap_length:
                path_info["default_gaps_num"] += 1
                path_info["default_gap_length_hpc"] += gap_length
            else:
                path_info["est_gaps_num"] += 1
                path_info["est_gap_length_hpc"] += gap_length
        else:
            contig, node_length = unitig_infos[member]
            contig_path.append(contig)
            path_members.add(member)
            path_info["nodes_num"] += 1
            path_info["total_path_length_hpc"] += node_length
            path_info["gapless_path_length_hpc"] += node_length
    
    if path_info["nodes_num"] > 1:
        path_info["path_type"] = "connected"

    path_info["contig_path"] = ",".join(contig_path)

    gapless_pct = path_info["gapless_path_length_hpc"] / path_info["total_path_length_hpc"]
    gapless_pct = round(gapless_pct * 100, 1)
    path_info["gapless_path_length_pct"] = gapless_pct

    unitig_records = []
    for node in path_members:
        node_info = dict(path_info)
        node_info["node"] = node
        unitig_records.append(node_info)

    return unitig_records


def main():

    args = parse_command_line()
    # create output folders
    for file_path in [args.out_table, args.out_table]:
        file_path.parent.mkdir(exist_ok=True, parents=True)

    unitigs = pd.read_csv(args.unitig_table, header=0, sep="\t")
    paths = pd.read_csv(args.rukki_paths, header=0, sep="\t")

    default_gap_length = args.default_gap
    node_infos = dict(
        (row.node, (row.contig, row.node_length_hpc))
        for row in unitigs.itertuples(index=False)
    )

    path_per_unitig = []
    for row in paths.itertuples(index=False):
        path_records = parse_rukki_path(
            default_gap_length,
            node_infos,
            row
        )
        path_per_unitig.extend(path_records)
    
    path_per_unitig = pd.DataFrame.from_records(path_per_unitig)
    unitigs = unitigs.merge(path_per_unitig, on="node", how="outer")
    unitigs["path_of_cc_length_pct"] = unitigs["gapless_path_length_hpc"] / unitigs["cc_length_hpc"]
    unitigs["path_of_cc_length_pct"] = (unitigs["path_of_cc_length_pct"] * 100).round(1)
    assert not pd.isnull(unitigs).any(axis=0).any()

    subset = unitigs.loc[unitigs["assigned"].isin(args.summary_subset), :].copy()
    subset.drop_duplicates(["path_id"], keep="first", inplace=True)
    subset.drop(
        [
            "node", "contig", "piece", "renamed",
            "node_length_hpc", "mat_marker",
            "pat_marker", "hap_group",
            "rukki_marker_sparsity",
            "ont_cov", "hifi_cov"
        ], axis=1, inplace=True
    )
    
    unitigs.to_csv(args.out_table, header=True, index=False, sep="\t")
    subset.to_csv(args.out_summary, header=True, index=False, sep="\t")

    return 0


if __name__ == "__main__":
    main()
