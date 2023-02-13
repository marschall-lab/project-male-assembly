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
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * 3}:59:00'
    run:
        import gzip as gz
        import collections as col
        import pandas as pd
        import numpy as np

        def karyo_loc(contig_name):   
            if "chrY" in contig_name:
                return "chrY"
            elif "chrX" in contig_name:
                return "chrX"
            else:
                return "auto"
        
        def compute_stats(contig_cov, contig_length):
            nzero_cov_bp = sum(contig_cov.values())
            zero_cov_bp = contig_length - nzero_cov_bp
            contig_cov[0] = zero_cov_bp
            contig_cov = pd.DataFrame.from_dict(
                contig_cov, orient="index", columns=["count"]
            )
            contig_cov.sort_index(inplace=True)
            cum_length = contig_cov.cumsum()
            median_cov = cum_length[cum_length > contig_length // 2].idxmin()
            mean_cov = round(
                np.average(contig_cov.index.tolist(), weights=contig_cov["count"].values),
                0
            )
            return int(mean_cov), int(median_cov)

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