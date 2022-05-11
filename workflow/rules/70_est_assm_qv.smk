
rule dump_short_read_kmers:
    """
    For now, this rule only supports
    a single input file
    """
    input:
        lambda wildcards: SAMPLE_DATA[wildcards.sample]['SHORT']
    output:
        'output/kmer_dump/{sample}.sr.yak'
    log:
        'log/output/kmer_dump/{sample}.sr.yak.log'
    benchmark:
        'rsrc/output/kmer_dump/{sample}.sr.yak.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, input, attempt: int(input.size_mb * 5 * attempt * 0.75),
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00'
    params:
        kmer_size = 31,
        bloom_size = 37,  # docs: this discards singletons
    shell:
        'yak count -k{params.kmer_size} -b{params.bloom_size} -t{threads} '
            '-o {output} {input} &> {log}'


rule estimate_assembly_qv:
    input:
        kmer_dump = 'output/kmer_dump/{sample}.sr.yak',
        fasta = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
    output:
        qv = 'output/eval/assembly_qv/{sample}.{hifi_type}.{ont_type}.na.wg.yak-qv.txt'
    log:
        'log/output/eval/assembly_qv/{sample}.{hifi_type}.{ont_type}.na.wg.yak-qv.log'
    benchmark:
        'rsrc/output/eval/assembly_qv/{sample}.{hifi_type}.{ont_type}.na.wg.yak-qv.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, input, attempt: int(input.size_mb * attempt * 0.75),
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00'
    params:
    shell:
        'yak qv -p -t{threads} {input.kmer_dump} {input.fasta} > {output.qv} 2> {log}'
