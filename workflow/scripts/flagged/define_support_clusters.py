#!/usr/bin/env python3

import argparse as argp
import collections as col
import hashlib as hl
import pathlib as pl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--regions",
        "-r",
        nargs="+",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="regions"
    )
    parser.add_argument(
        "--depths",
        "-d",
        nargs="+",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="depths"
    )
    parser.add_argument(
        "--merged",
        "-m",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="merged"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        required=True
    )
    args = parser.parse_args()
    return args


def add_coverage(flanked_regions, depth_table, read_label):
    
    rd = pd.read_csv(depth_table, sep="\t", header=None, names=["contig", "pos", "depth"])
    # NB: pos is 1-based > make zero-based
    rd["pos"] -= 1
    depths = []
    for contig, regions in flanked_regions.groupby("contig"):
        sub_ctg = rd.loc[rd["contig"] == contig, :]
        for region in regions.itertuples():
            region_length = region.end - region.start
            select_start = region.start <= sub_ctg["pos"]
            select_end = region.end > sub_ctg["pos"]
            cov_values = sub_ctg.loc[select_start & select_end, "depth"]
            if cov_values.size != region_length:
                # boundary effect... make sure it's flanking
                assert region.region_type == 'flank_3p'
                cov_values = sub_ctg.loc[select_start, "depth"]
            depths.append(
                (
                    region.name,
                    round(cov_values.mean(), 1),
                    round(cov_values.median(), 1)
                )
            )
    depths = pd.DataFrame.from_records(
        depths, columns=["name", f"{read_label}_mean_cov", f"{read_label}_median_cov"]
    )
    flanked_regions = flanked_regions.merge(depths, on="name")
    assert not pd.isnull(flanked_regions).any(axis=1).any()
    return flanked_regions


def combine_all_regions(region_files, depth_files):

    read_labels = {"HIFIRW": "hifi", "ONTUL": "ont"}

    all_regions = []
    for region_file in region_files:
        region_source = region_file.name.split(".")[-2]
        regions = pd.read_csv(region_file, sep="\t", header=0)   
        matched_depths = [fp for fp in depth_files if region_source in fp.name]
        # note that the above is empty list for SNV regions
        for depth_file in matched_depths:
            read_type = depth_file.name.split("_aln-to_")[0]
            read_type = read_type.split(".")[1]
            read_label = read_labels[read_type]
            regions = add_coverage(regions, depth_file, read_label)
        all_regions.append(regions)
    
    all_regions = pd.concat(all_regions, axis=0, ignore_index=False)
    all_regions.fillna(-1, inplace=True)
    all_regions.sort_values(["contig", "start", "end"], inplace=True)

    assert not pd.isnull(all_regions).any(axis=1).any()

    return all_regions


def identify_support_clusters(regions, merged):

    tool_label_map = {
        "NucFreq": "num_nucfreq_regions",
        "VerityMap": "num_veritymap_regions",
        "DeepVariant": "num_het_snv_hifi",
        "PEPPER": "num_het_snv_ont",
    }
    assert all(tool in tool_label_map for tool in regions["software"].unique())

    cluster_stats = []
    for cluster in merged.itertuples(index=False):
        region_ids = sorted(cluster.region_ids.split(","))
        if len(region_ids) == 1:
            # nothing overlaps
            continue
        origins = regions.loc[regions["name"].isin(region_ids), :]
        if origins["software"].nunique() == 1:
            # everything merged came from the same tool,
            # not interesting; coverage will show if region is problematic
            continue
        if origins["software"].isin(["DeepVariant", "NucFreq", "PEPPER"]).all():
            # NucFreq and variant callers essentially look
            # for the same thing, so not interesting
            continue
        tool_counts = origins["software"].value_counts()
        tool_support = []
        total_snvs = 0
        for tool in tool_label_map.keys():
            try:
                tool_count = tool_counts[tool]
            except KeyError:
                tool_count = 0
            tool_support.append(tool_count)
            if tool in ["DeepVariant", "PEPPER"]:
                total_snvs += tool_count
        
        start = origins["start"].min()
        end = origins["end"].max()
        cluster_id = hl.md5("".join(region_ids).encode("utf-8")).hexdigest()
        cluster_length = end - start
        snv_density = total_snvs / (cluster_length / 1000)
        for region in region_ids:
            record = [region, cluster_id, start, end, cluster_length, snv_density]
            record.extend(tool_support)
            cluster_stats.append(record)

    stats_columns = [
        "name", "cluster_id", "cluster_start", "cluster_end", "cluster_span", "snv_density_kbp"
    ]
    stats_columns.extend(list(tool_label_map.values()))
    cluster_stats = pd.DataFrame.from_records(
        cluster_stats, 
        columns=stats_columns
    )

    assert not pd.isnull(cluster_stats).any(axis=1).any()

    return cluster_stats


def main():

    args = parse_command_line()

    # contains everything, including flanking regions
    regions = combine_all_regions(args.regions, args.depths)
    
    merged_regions = pd.read_csv(
        args.merged,
        sep="\t",
        header=None,
        names=["contig", "start", "end", "region_ids"]
    )

    # called by tools / not flanking
    origin_regions = regions.loc[regions["region_type"] == "origin", :]

    support_clusters = identify_support_clusters(origin_regions, merged_regions)

    regions = regions.merge(support_clusters, on="name", how="outer")
    regions["cluster_id"].fillna("singleton", inplace=True)
    for column in support_clusters.columns:
        if column in ["name", "cluster_id"]:
            continue
        regions[column].fillna(-1, inplace=True)
        if column in ["cluster_start", "cluster_end", "cluster_span"]:
            regions[column] = regions[column].astype(int)
        
    assert not pd.isnull(regions).any(axis=1).any()

    regions.sort_values(["contig", "start", "end"], inplace=True)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    regions.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
