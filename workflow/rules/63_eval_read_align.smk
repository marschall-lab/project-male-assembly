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


rule run_all_read_depth:
    input:
        depths = expand(
            'output/eval/read_cov/wg/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.minmq{minmapq}.depth.tsv.gz',
            sample=COMPLETE_SAMPLES,
            other_reads=["HIFIRW", "ONTUL"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            minmapq=[0, 1, 10]
        )