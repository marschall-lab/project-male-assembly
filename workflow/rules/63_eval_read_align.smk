"""
Module added for revision to
produce different views on
genome-wide read coverage
(for de-novo assemblies)
"""


rule dump_read_to_assm_coverage:
    """
    NB: samtools depth skips reads ...

        > By default, reads that have any of the flags
        > UNMAP, SECONDARY, QCFAIL, or DUP set are skipped.
        (see man pages)

    """
    input:
        bed = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.{genome}.ctg-500kbp.bed',
        bam = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai'
    output:
        depth = 'output/eval/read_cov/wg/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{genome}.minmq{minmapq}.depth.tsv.gz'
    wildcard_constraints:
        genome = "(noYHET|wg)"
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt}:59:00'
    shell:
        'samtools depth --min-MQ {wildcards.minmapq} -l 5000 '
        '-b {input.bed} {input.bam} | pigz -p 2 > {output.depth}'


rule agg_read_to_assm_coverage:
    input:
        contigs = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.bed',
        no_hets = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.noYHET.ctg-500kbp.bed',
        depth = 'output/eval/read_cov/wg/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{genome}.minmq{minmapq}.depth.tsv.gz',
    output:
        table = 'output/eval/read_cov/stats/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{genome}.minmq{minmapq}.stats.tsv',
    wildcard_constraints:
        genome = "(noYHET|wg)"
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * 3}:59:00'
    run:
        import gzip as gz
        import collections as col
        import pandas as pd
        import numpy as np

        NO_HET = wildcards.genome == "noYHET"

        def karyo_loc(contig_name):   
            if "chrY" in contig_name:
                return "chrY"
            elif "chrX" in contig_name:
                return "chrX"
            else:
                return "auto"
        
        def compute_stats(contig_cov, contig_length):
            nzero_cov_bp = sum(contig_cov.values())
            assert nzero_cov_bp >= 0
            assert nzero_cov_bp <= contig_length, f"NZ-cov {nzero_cov_bp} / CL {contig_length}"
            zero_cov_bp = contig_length - nzero_cov_bp
            assert zero_cov_bp >= 0
            contig_cov[0] += zero_cov_bp
            assert sum(contig_cov.values()) == contig_length
            contig_cov = pd.DataFrame.from_dict(
                contig_cov, orient="index", columns=["count"]
            )
            contig_cov.sort_index(inplace=True)
            cum_length = contig_cov.cumsum()
            median_cov = cum_length[cum_length > contig_length // 2].idxmin().values[0]
            mean_cov = round(
                np.average(contig_cov.index.tolist(), weights=contig_cov["count"].values),
                0
            )
            if pd.isnull(mean_cov) or pd.isnull(median_cov):
                print("=====ERROR=====")
                print(contig_cov)
                print("=====ERROR=====")
                raise ValueError(
                    f"Mean or median undefined: CL {contig_length} / AVG {mean_cov} / MED {median_cov}"
                )
            return int(mean_cov), int(median_cov)

        if NO_HET:
            contigs = pd.read_csv(
                input.no_hets, sep="\t",
                names=["contig", "start", "end"],
            )
            contigs["length"] = contigs["end"] - contigs["start"]
            contigs = contigs.groupby("contig")["length"].sum()
            contigs = contigs.reset_index(inplace=False, drop=False)
        else:
            contigs = pd.read_csv(
                input.contigs, sep="\t",
                names=["contig", "start", "length"],
                usecols=["contig", "length"]
            )
        contigs["location"] = contigs["contig"].apply(karyo_loc)
        total_size = contigs["length"].sum()
        if NO_HET:
            # NB: when removing all HET regions from chrY,
            # the contig size filter > 500 kbp was already
            # applied
            proc_size = total_size
        else:
            proc_size = contigs.loc[contigs["length"] >= 500000, "length"].sum()
        karyo_size_total = contigs.groupby("location")["length"].sum()
        if NO_HET:
            # same as above
            karyo_size_proc = karyo_size_total
        else:
            karyo_size_proc = contigs.loc[contigs["length"] >= 500000, :].groupby("location")["length"].sum()
        proc_pct = round(proc_size / total_size * 100, 2)

        if NO_HET:
            cov_stats = []
            stats_context = "no_het"
        else:
            stats_context = "global"
            cov_stats = [
                ("wg", "size_bp", "global", total_size),
                ("wg", "proc_bp", "global", proc_size),
                ("wg", "proc_pct", "global", proc_pct),   
            ]
            

        for loc, size_bp in karyo_size_total.items():
            cov_stats.append((loc, "size_bp", stats_context, size_bp))
            proc_size = karyo_size_proc.at[loc]
            proc_pct = round(proc_size / size_bp * 100, 2)
            cov_stats.append((loc, "proc_bp", stats_context, proc_size))
            cov_stats.append((loc, "proc_pct", stats_context, proc_pct))

        contig_cov = col.Counter()
        current_contig = None
        with gz.open(input.depth, "rt") as table:
            for line in table:
                ctg, pos, cov = line.split()
                if ctg != current_contig and current_contig is not None:
                    chrom_loc, contig_length = contigs.loc[
                        contigs["contig"] == current_contig, ["location", "length"]
                    ].values[0]
                    mean_cov, median_cov = compute_stats(contig_cov, contig_length)
                    cov_stats.append((chrom_loc, "mean_cov", current_contig, mean_cov))
                    cov_stats.append((chrom_loc, "median_cov", current_contig, median_cov))
                    # reset structs
                    contig_cov = col.Counter()
                    contig_length = 0
                current_contig = ctg
                contig_cov[int(cov)] += 1

        chrom_loc, contig_length = contigs.loc[
            contigs["contig"] == current_contig, ["location", "length"]
            ].values[0]
        mean_cov, median_cov = compute_stats(contig_cov, contig_length)
        cov_stats.append((chrom_loc, "mean_cov", current_contig, mean_cov))
        cov_stats.append((chrom_loc, "median_cov", current_contig, median_cov))          
            
        cov_stats = pd.DataFrame.from_records(
            cov_stats,
            columns=["location", "statistic", "context", "value"]
        )
        cov_stats.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule combine_read_depth_stats:
    input:
        tables = expand(
            'output/eval/read_cov/stats/{sample}.{other_reads}_aln-to_{{hifi_type}}.{{ont_type}}.na.{genome}.minmq{minmapq}.stats.tsv',
            sample=COMPLETE_SAMPLES,
            other_reads=["HIFIRW", "ONTUL"],
            genome=["noYHET", "wg"],
            minmapq=[0, 10]
        )
    output:
        table = 'output/eval/read_cov/stats/ALL-SAMPLES.READS_aln-to_{hifi_type}.{ont_type}.na.wg.cov-stats.tsv'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import pathlib as pl

        merged = []
        for table_file in input.tables:
            prefix, suffix = pl.Path(table_file).name.split("_aln-to_")
            sample, readset = prefix.split(".")
            mapq = int(suffix.split(".")[-3].strip("minmq"))
            input_set = suffix.split(".")[-4]
            assert input_set in ["wg", "noYHET"]
            table = pd.read_csv(table_file, sep="\t", header=0)
            table["sample"] = sample
            table["reads"] = readset
            table["min_mapq"] = mapq
            table["input_set"] = input_set
            merged.append(table)
        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["sample", "context", "input_set", "reads", "min_mapq"])
        merged.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_read_depth:
    input:
        merged = expand(
            'output/eval/read_cov/stats/ALL-SAMPLES.READS_aln-to_{hifi_type}.{ont_type}.na.wg.cov-stats.tsv',
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
        )