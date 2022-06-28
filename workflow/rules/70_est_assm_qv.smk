
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
        mem_mb = lambda wildcards, attempt: 65536 + 32768 * attempt,
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
        mem_mb = lambda wildcards, input, attempt: int(input.size_mb * attempt * 0.8),
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00'
    params:
    shell:
        'yak qv -p -t{threads} {input.kmer_dump} {input.fasta} > {output.qv} 2> {log}'


# add alternative method
# run Merqury on whole-genome assembly

rule dump_meryl_kmer_db:
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample][wildcards.reads]
    output:
        db = directory('output/kmer_dump/{sample}.{reads}.k{kmer}.meryl')
    log:
        'log/output/kmer_dump/{sample}.{reads}.k{kmer}.meryl.log'
    benchmark:
        'rsrc/output/kmer_dump/{sample}.{reads}.k{kmer}.meryl.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, input, attempt: int(input.size_mb * 2 * attempt),
        mem_gb = lambda wildcards, input, attempt: int(input.size_mb * 2 * attempt / 1024),
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00'
    params:
        kmer_size = lambda wildcards: int(wildcards.kmer)
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_gb} output {output.db} {input.reads} &> {log}'


rule run_merqury:
    """
    Merqury was added as an alternative to yak,
    but because of the mess it creates with its
    input symlinking strategy, it is essentially
    not usable. Generating the output of this rule
    is currently not triggered when running the
    pipeline with targets "run_all" or
    "populate_share".
    """
    input:
        meryl_db = 'output/kmer_dump/{sample}.{other_reads}.k{kmer}.meryl',
        wg_assm = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
    output:
        chk = 'output/eval/merqury/{sample}.{hifi_type}.{ont_type}.na.wg_{other_reads}.k{kmer}.ok'
    log:
        'log/output/eval/merqury/{sample}.{hifi_type}.{ont_type}.na.wg_{other_reads}.k{kmer}.merqury.log'
    benchmark:
        'rsrc/output/eval/merqury/{sample}.{hifi_type}.{ont_type}.na.wg_{other_reads}.k{kmer}.merqury.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, input, attempt: int(input.size_mb * 2 * attempt),
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00'
    params:
        outdir = lambda wildcards, output: output.chk.rsplit('.', 1)[-2]
    shell:
        '$MERQURY/merqury.sh {input.meryl_db} {input.wg_assm} {params.outdir} &> {log}'
            ' && '
        'touch {output.chk}'
