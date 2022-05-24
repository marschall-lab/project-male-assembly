
rule merge_read_sets:
    """
    VerityMap only supports a single
    input file for the reads...
    """
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample][wildcards.read_type]
    output:
        reads = temp('temp/merged_reads/{sample}_{read_type}.fa.gz')
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        walltime = lambda wildcards, attempt: f'{23*attempt}:59:00'
    shell:
        'seqtk seq -A -C {input.reads} | pigz -p {threads} --best > {output.reads}'


rule run_veritymap:
    """
    Jobs failed apparently after restart b/c of a wrongly
    resolved path from Singularity POV. Presumably, when
    trying to resuse existing results, VerityMap performs
    some copy operations to that end, which then fail.
    Add "no-reuse" statement and empty output dir
    """
    input:
        reads = 'temp/merged_reads/{sample}_{read_type}.fa.gz',
        assembly = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
    output:
        chk = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.chk'
    log:
        'log/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.log'
    benchmark:
        'rsrc/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.rsrc'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap"]}'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 65536 * attempt,
        walltime = lambda wildcards, attempt: f'{23*attempt}:59:00'
    params:
        preset = lambda wildcards: {'HIFIRW': 'hifi', 'ONTUL': 'ont'}[wildcards.read_type],
        outdir = lambda wildcards, output: output.chk.rsplit('.', 2)[0]
    shell:
        'rm -rf {params.outdir} && '
        'python /repos/VerityMap/veritymap/main.py --no-reuse --reads {input.reads} -t {threads} '
        '-d {params.preset} -l {wildcards.sample} -o {params.outdir} {input.assembly} &> {log}'
            ' && '
        'touch {output.chk}'
