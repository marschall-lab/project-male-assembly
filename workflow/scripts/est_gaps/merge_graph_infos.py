#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl

import pandas as pd
import networkx as nx


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--graph",
        "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="graph",
        help="Path to Verkko output graph (GFA) file."
    )

    parser.add_argument(
        "--scf-layout",
        "-l",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="layout",
        help="Path to Verkko layout file ('unitig-popped.layout.scfmap')."
    )

    parser.add_argument(
        "--contig-names",
        "-n",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs=2,
        dest="name_files",
        help="Path to contig (re-) name files (for chrX and chrY)."
    )

    parser.add_argument(
        "--hifi-node-cov",
        "-hc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="hifi_cov",
        help="Path to node HIFI coverage file."
    )

    parser.add_argument(
        "--ont-node-cov",
        "-oc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="ont_cov",
        help="Path to node ONT coverage file."
    )

    parser.add_argument(
        "--out-unitigs",
        "-o",
        type=lambda x: pl.Path(x).resolve(),
        dest="unitigs",
        help="Path to output table listing all unitigs with merged infos (TSV)."
    )

    parser.add_argument(
        "--out-coloring",
        "-c",
        type=lambda x: pl.Path(x).resolve(),
        dest="coloring",
        help="Path to output coloring for Rukki (TSV)."
    )

    parser.add_argument(
        "--out-conncomp",
        "-p",
        type=lambda x: pl.Path(x).resolve(),
        dest="components",
        help="Path to output connected components (TSV)."
    )

    args = parser.parse_args()

    return args


def read_name_file(file_path, contig_to_chrom):
    """Works via side effect; directly manipulates
    contig to chrom mapping dict
    """

    with open(file_path, "r") as table:
        for line in table:
            full_name = line.strip()
            components = full_name.split(".")
            chrom = components[0]
            assert chrom in ["chrX", "chrY"]
            contig_name = components[-2]
            assert "unassigned" in contig_name
            sample = components[-1]
            contig_to_chrom[contig_name] = {
                "chrom": chrom,
                "renamed": full_name,
                "sample": sample
            }
    return None


def read_coverage_file(file_path, data_source):

    cov_table = pd.read_csv(
        file_path,
        header=0,
        sep="\t"
    )
    cov_table[f"{data_source}_cov"] = cov_table["coverage"].round(2)
    cov_table.drop("coverage", axis=1, inplace=True)
    return cov_table


def read_layout_file(file_path):

    current_utg = None
    current_ctg = None
    tig_mapping = dict()
    with open(file_path, "r") as listing:
        for ln, line in enumerate(listing, start=1):
            if line.startswith("path"):
                _, contig, unitig = line.strip().split()
                current_ctg = contig
                current_utg = unitig
            elif line.startswith("piece"):
                piece = line.strip()
                assert current_utg not in tig_mapping
                assert current_utg is not None
                tig_mapping[current_utg] = {
                    "contig": current_ctg,
                    "piece": piece
                }
            elif line.startswith("end"):
                current_utg = None
                current_ctg = None
            else:
                raise ValueError(f"{ln}: {line.strip()}")
    return tig_mapping


def build_graph(graph_file):

    graph = nx.Graph()
    node_lengths = col.Counter()
    with open(graph_file, "r") as gfa:
        for ln, line in enumerate(gfa, start=1):
            if line.startswith("S"):
                columns = line.strip().split()
                node = columns[1]
                node_length = len(columns[2])
                node_lengths[node] = node_length
                graph.add_node(node)
            elif line.startswith("L"):
                columns = line.strip().split()
                node_a = columns[1]
                node_b = columns[3]
                graph.add_edge(node_a, node_b)
            else:
                raise ValueError(f"{ln}: {line.strip()}")
    return graph, node_lengths


def process_connected_components(graph, node_lengths, layout, ctg_names):

    unitig_records = []
    cc_records = []
    for cc_id, cc_nodes in enumerate(nx.connected_components(graph), start=1):
        cc = graph.subgraph(cc_nodes).copy()
        size_of_cc = len(cc.nodes)
        length_of_cc = sum(node_lengths[n] for n in cc.nodes)
        cc_record = {
            "cc_id": f"CC{cc_id}",
            "cc_size": size_of_cc,
            "cc_length_hpc": length_of_cc
        }
        cc_records.append(cc_record)
        for node in cc.nodes:
            # node = unitig
            unitig_record = dict(cc_record)
            unitig_record["node"] = node
            node_length = node_lengths[node]
            unitig_record["node_length_hpc"] = node_length
            marker_count = max(100, int(node_length/1000))
            # 10000 is the default in rukki
            marker_sparsity = int(round(node_length/10000, 0))
            try:
                layout_info = layout[node]
                contig = layout_info["contig"]
                piece = layout_info["piece"]
            except KeyError:
                contig = "no-layout"
                piece = "no-piece"
            unitig_record["contig"] = contig
            unitig_record["piece"] = piece
            try:
                contig_info = ctg_names[contig]
                chrom = contig_info["chrom"]
                renamed = contig_info["renamed"]
                sample = contig_info["sample"]
                mat_marker = 0
                pat_marker = 0
                hap_group = 0
                if chrom == "chrX":
                    mat_marker = marker_count
                    hap_group = 1
                if chrom == "chrY":
                    pat_marker = marker_count
                    hap_group = 2
            except KeyError:
                chrom = "chrUN"
                renamed = "UNK"
                sample = "UNK"
                mat_marker = 0
                pat_marker = 0
                hap_group = 0
            unitig_record["chrom"] = chrom
            unitig_record["mat_marker"] = mat_marker
            unitig_record["pat_marker"] = pat_marker
            unitig_record["hap_group"] = hap_group
            unitig_record["rukki_marker_sparsity"] = marker_sparsity
            unitig_record["renamed"] = renamed
            unitig_record["sample"] = sample
            unitig_records.append(unitig_record)

    cc_df = pd.DataFrame.from_records(cc_records)
    assert not pd.isnull(cc_df).all(axis=0).all()
    utg_df = pd.DataFrame.from_records(unitig_records)
    assert not pd.isnull(utg_df).all(axis=0).all()

    return cc_df, utg_df



def main():

    args = parse_command_line()
    # create output folders
    for file_path in [args.unitigs, args.coloring, args.components]:
        file_path.parent.mkdir(exist_ok=True, parents=True)
    graph, node_lengths = build_graph(args.graph)
    # create tig mapping
    layout = read_layout_file(args.layout)
    contig_to_chrom = dict()
    [read_name_file(fp, contig_to_chrom) for fp in args.name_files]

    conncomp, unitigs = process_connected_components(
        graph, node_lengths, layout, contig_to_chrom
    )
    hifi_cov = read_coverage_file(args.hifi_cov, "hifi")
    ont_cov = read_coverage_file(args.ont_cov, "ont")
    unitigs = unitigs.merge(hifi_cov, on="node", how="outer")
    unitigs["hifi_cov"].fillna(0, inplace=True)

    unitigs = unitigs.merge(ont_cov, on="node", how="outer")
    unitigs["ont_cov"].fillna(0, inplace=True)
    assert not pd.isnull(unitigs).all(axis=0).all()
    unitigs.sort_values(
        ["hap_group", "cc_length_hpc", "node_length_hpc"],
        ascending=[False, False, False],
        inplace=True
    )

    conncomp.to_csv(args.components, sep="\t", header=True, index=False)
    unitigs.to_csv(args.unitigs, sep="\t", header=True, index=False)
    unitigs[["node", "mat_marker", "pat_marker"]].to_csv(
        args.coloring, sep="\t", header=True, index=False
    )


if __name__ == "__main__":
    main()
