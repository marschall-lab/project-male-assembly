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
        bed = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.ctg-500kbp.bed',
        bam = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai'
    output:
        depth = 'output/eval/read_cov/wg/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.minmq{minmapq}.depth.tsv.gz'
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
        depth = 'output/eval/read_cov/wg/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.minmq{minmapq}.depth.tsv.gz',
    output:
        table = 'output/eval/read_cov/stats/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.minmq{minmapq}.stats.tsv',
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt}:59:00'
    run:
        import pandas as pd
        import numpy as np

        def karyo_loc(contig_name):   
            if "chrY" in contig_name:
                return "chrY"
            elif "chrX" in contig_name:
                return "chrX"
            else:
                return "auto"

        contigs = pd.read_csv(
            input.contigs, sep="\t",
            names=["contig", "start", "length"],
            usecols=["contig", "length"]
        )
        contigs["location"] = contigs["contig"].apply(karyo_loc)
        total_size = contigs["length"].sum()
        proc_size = contigs.loc[contigs["length"] >= 500000, "length"].sum()
        karyo_size_total = contigs.groupby("location")["length"].sum()
        karyo_size_proc = contigs.loc[contigs["length"] >= 500000, :].groupby("location")["length"].sum()
        proc_pct = round(proc_size / total_size * 100, 2)

        cov_stats = [
            ("wg", "size_bp", "global", total_size),
            ("wg", "proc_bp", "global", proc_size),
            ("wg", "proc_pct", "global", proc_pct),   
        ]

        for loc, size_bp in karyo_size_total.items():
            cov_stats.append((loc, "size_bp", "global", size_bp))
            proc_size = karyo_size_proc.at[loc]
            proc_pct = round(proc_size / size_bp * 100, 2)
            cov_stats.append((loc, "proc_bp", "global", proc_size))
            cov_stats.append((loc, "proc_pct", "global", proc_pct))
 
        depth_table = pd.read_csv(
            input.depth,
            sep="\t",
            names=["contig", "pos", "cov"],
            usecols=["contig", "cov"]
        )

        for contig, coverage in depth_table.groupby("contig"):
            chrom_loc, contig_length = contigs.loc[contigs["contig"] == contig, ["location", "length"]].values[0]
            covered_positions = coverage.shape[0]
            cov_zero = contig_length - covered_positions
            cov_counts = coverage["cov"].value_counts()
            # NB: this 0 refers to the index of the Series
            if 0 in cov_counts:
                # NB: this is at index position 0,
                # not at first position of Series
                cov_counts[0] += cov_zero
            else:
                cov_zero = pd.Series([cov_zero], index=[0])
                cov_counts = pd.concat([cov_counts, cov_zero], ignore_index=False)
            cov_counts.sort_index(inplace=True)
            cum_length = cov_counts.cumsum()
            median_cov = cum_length[cum_length > contig_length // 2].idxmin()
            mean_cov = np.average(cov_counts.index.tolist(), weights=cov_counts.values)
            cov_stats.append((chrom_loc, "mean_cov", contig, mean_cov))
            cov_stats.append((chrom_loc, "median_cov", contig, median_cov))
            
        cov_stats = pd.DataFrame.from_records(cov_stats, columns=["location", "statistic", "context", "value"])



rule run_all_read_depth:
    input:
        stats = expand(
            'output/eval/read_cov/stats/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.minmq{minmapq}.stats.tsv',
            sample=COMPLETE_SAMPLES,
            other_reads=["HIFIRW", "ONTUL"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            minmapq=[0, 1, 10]
        )