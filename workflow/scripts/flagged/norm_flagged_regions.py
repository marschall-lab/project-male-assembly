#!/usr/bin/env python

import argparse as argp
import collections as col
import csv
import dataclasses as dcl
import functools as fnt
import hashlib as hl
import pathlib as pl

import pandas as pd

def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--table",
        "-t",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="table",
        help="Path to input table"
    )
    parser.add_argument(
        "--curated",
        "-c",
        type=lambda x: pl.Path(x).resolve(strict=True),
        default=None,
        dest="curated",
        help="Path to expert curated table."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=lambda x: pl.Path(x).resolve(strict=False)
    )
    args = parser.parse_args()
    return args


@dcl.dataclass
class GenomicRegion:
    contig: str
    start: int
    end: int
    name: str
    sample: str
    est_size: int
    read_source: str
    software: str
    region_type: str
    region_source: str
    region_start: int
    region_end: int
    curated: int = 0
    hifi_support: int = 0
    ont_support: int = 0
    _key_order: tuple = (
        "contig", "start", "end", "name",
        "est_size", "sample", "region_type",
        "region_source", "curated", "ont_support",
        "hifi_support", "software", "read_source",
        "region_start", "region_end"
    )
    _add_fields: tuple = None

    def __init__(self, 
        contig, start, end, sample, est_size,
        read_source, software,
        region_type="origin", region_source=None,
        add_fields=None):
        
        self.contig = contig
        self.region_start = int(start)
        self.region_end = int(end)
        assert self.region_start >= 0, f"{contig} - {start} - {end}"
        assert self.region_start <= self.region_end, f"{contig} - {start} - {end}"
        self.est_size = int(est_size)
        self._set_region_boundary()
        self._set_name()
        self.sample = sample
        self.read_source = read_source
        self.software = software
        self.region_type = region_type
        if region_source is None:
            assert self.region_type == "origin"
            self.region_source = self.name
        else:
            self.region_source = region_source
        self._add_fields = add_fields
        return None
   
    def _set_name(self):
        
        parts = ".".join([self.contig, str(self.start), str(self.end)])
        name = hl.md5(parts.encode("utf-8")).hexdigest()
        self.name = name
        return None
    
    def _set_region_boundary(self):
        
        region_size = abs(self.est_size)
        if self.region_end - self.region_start < 2:
            self.start = max(0, self.region_start - region_size // 2)
            self.end = self.region_end + region_size // 2
        else:
            self.start = self.region_start
            self.end = self.region_end
        assert self.start >= 0
        assert self.start < self.end
        return None
    
    def make_flanking(self, side):
        
        length = self.end - self.start
        if side in ["5p", "up", "left", "before"]:
            new_start = max(0, self.start - length)
            new_end = self.start
            region_type = "flank_5p"
        elif side in ["3p", "down", "right", "after"]:
            new_start = max(0, self.end)
            new_end = self.end + length
            region_type = "flank_3p"
        else:
            raise ValueError(side)
        assert new_start < new_end
        new_region = GenomicRegion(
            self.contig, new_start, new_end, self.sample,
            length, self.read_source, self.software,
            region_type=region_type, region_source=self.name
            
        )
        return new_region

    def check_curated(self, annotation):
        
        try:
            expert_label = annotation[self.name]
        except KeyError:
            pass
        else:
            self.curated = 1
            self.ont_support = expert_label["ont"]
            self.hifi_support = expert_label["hifi"]
        return None
    
    def to_dict(self):
        
        gr_dict = col.OrderedDict()
        for key in self._key_order:
            gr_dict[key] = getattr(self, key)
        if self._add_fields is not None:
            for key, value in self._add_fields:
                gr_dict[key] = value
        return gr_dict
    
    @classmethod
    def get_header(cls):
        return cls._key_order


def parse_norm_row(sample, row):

    assert sample == row["sample"]
    gr = GenomicRegion(
        row["chrom"], row["start"], row["end"],
        sample, row["est_size"], row["reads"],
        row["source"]
    )
    return gr


def parse_nucfreq_row(sample, row):

    software = "NucFreq"
    start = int(row["start"])
    end = int(row["end"])
    est_size = end - start
    add_fields = (
        ("nucfreq_median_het_ratio", row["median_het_ratio"]),
        ("nucfreq_num_hets", row["num_hets"])
    )
    gr = GenomicRegion(
        row["contig"], start, end,
        sample, est_size, "HIFIRW",
        software, add_fields=add_fields
    )
    return gr


def _relabel(label):
    if label == 1:
        return 1
    if label == 0:
        return -1
    if label == -1:
        return 0
    raise label
        
def _make_name(row):
    parts = f"{row.contig}.{row.start}.{row.end}"
    name = hl.md5(parts.encode("utf-8")).hexdigest()
    return name


def process_expert_label_file(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)
    df.drop_duplicates(inplace=True)
    df["label_hifi"] = df["hifi_support"].apply(_relabel)
    df["label_ont"] = df["ont_support"].apply(_relabel)
    df["name"] = df.apply(_make_name, axis=1)
    lut = dict(
        (row.name, {"hifi": row.label_hifi, "ont": row.label_ont})
        for row in df.itertuples(index=False)
    )
    assert len(lut) == df.shape[0]
    return lut


def main():

    args = parse_command_line()

    file_sample = args.table.name.split(".")[0]
    if args.table.name.endswith("het-regions.tsv"):
        row_parser = fnt.partial(parse_nucfreq_row, file_sample)
    else:
        row_parser = fnt.partial(parse_norm_row, file_sample)

    if "pr-errors.tsv" in args.table.name or "dv-errors.tsv" in args.table.name:
        flank_region = False
    else:
        flank_region = True

    expert_labels = None
    if args.curated is not None:
        expert_labels = process_expert_label_file(args.curated)

    table_regions = []
    with open(args.table, "r", newline="") as table:
        header = table.readline().strip().split()
        rows = csv.DictReader(table, fieldnames=header, delimiter="\t")
        for row in rows:
            region = row_parser(row)
            if expert_labels is not None:
                region.check_curated(expert_labels)
            table_regions.append(region)
            if flank_region:
                left = region.make_flanking("left")
                table_regions.append(left)
                right = region.make_flanking("right")
                table_regions.append(right)

    df = pd.DataFrame.from_records(
        [r.to_dict() for r in table_regions]
    )
    if df.empty:
        # this happens for high-quality samples,
        # and commonly b/c does not call variants
        header = GenomicRegion.get_header()
        df = pd.DataFrame(columns=header)
    df.sort_values(["contig", "start", "end"], inplace=True)
    df.fillna(0, inplace=True)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(args.output, header=True, index=False, sep="\t")
    bed_output = args.output.with_suffix(".bed")
    with open(bed_output, "w") as dump:
        _ = dump.write("#")
        df[["contig", "start", "end", "name"]].to_csv(
            dump, header=True, index=False, sep="\t"
        )

    return 0
        



if __name__ == "__main__":
    main()
