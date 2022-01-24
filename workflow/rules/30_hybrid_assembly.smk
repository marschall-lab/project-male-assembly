
rule run_verkko_targeted_assembly:
    input:
        ont = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{ont_type}.{chrom}-reads.{mapq}.fasta.gz',
        hifi = select_hifi_reads
    output:
        wd = directory('output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}'),
        complete = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.ok'
    log:
        'log/output/hybrid/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.log'
    benchmark:
        'rsrc/output/hybrid/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.rsrc'
    conda:
        '../envs/verkko.yaml'
    threads: config['num_cpu_high']
    resources:
        mem_mb = lambda wildcards, attempt: 110592 * attempt,
        walltime = lambda wildcards, attempt: f'{35 * attempt}:59:00'
    params:
        run_correction = lambda wildcards: '--no-correction' if wildcards.hifi_type in ['HIFIEC', 'ONTEC', 'OHEC'] else '',
        high_maxk = lambda wildcards: '--max-k 100000' if wildcards.hifi_type in ['ONTEC', 'OHEC'] else ''
    shell:
        'verkko --local --hifi {input.hifi} --nano {input.ont} -d {output.wd} '
            '{params.run_correction} {params.high_maxk} --threads {threads} &> {log} '
            ' && touch {output.complete}'
