
localrules: aggregate_yak_assembly_qv

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
    shell:
        'yak qv -p -t{threads} {input.kmer_dump} {input.fasta} > {output.qv} 2> {log}'


rule aggregate_yak_assembly_qv:
    input:
        tables = expand(
            'output/eval/assembly_qv/{sample}.{{hifi_type}}.{{ont_type}}.na.wg.yak-qv.txt',
            sample=COMPLETE_SR_SAMPLES
        )
    output:
        table = 'output/eval/assembly_qv/SAMPLES.{hifi_type}.{ont_type}.na.wg.yak-qv.tsv',
    run:
        import pathlib as pl

        assembly_qv = []
        for file_path in input.tables:
            sample_name = pl.Path(file_path).stem.split('.')[0]
            is_qc_only = 1 if sample_name in QC_SAMPLES else 0
            with open(file_path, 'r') as table:
                for line in table:
                    if not line.startswith('QV'):
                        continue
                    _, raw, adjusted = line.strip().split()
                    assembly_qv.append((sample_name, 'wg', 'yak', adjusted, str(is_qc_only)))

        with open(output.table, 'w') as table:
            _ = table.write('sample\tassembly\test_method\tqv\tqc_only_assembly\n')
            for record in sorted(assembly_qv):
                _ = table.write('\t'.join(record) + '\n')

    # END OF RUN BLOCK


# attempt to get a reasonable QV estimate using short-reads just for chrY

rule map_short_reads_augmented_assembly:
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample]["SHORT"],
        assm_idx = multiext(
            "output/output/eval/chry_qv/aug_assm/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.T2TY.idx/{sample}",
            ".64.amb", ".64.ann", ".64.bwt", ".64.pac", ".64.sa"
        )
    output:
        bam = "output/output/eval/chry_qv/aln_assm_sr/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.T2TY.short.bam",
        bai = "output/output/eval/chry_qv/aln_assm_sr/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.T2TY.short.bam.bai",
    log:
        "log/output/output/eval/chry_qv/aln_assm_sr/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.T2TY.short.bwa.log"
    benchmark:
        "rsrc/output/output/eval/chry_qv/aln_assm_sr/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.T2TY.short.bwa.rsrc"
    conda:
        "../envs/biotools.yaml"
    params:
        prefix = lambda wildcards, output: output.idx[0].rsplit(".", 2)[0]
    threads: config["num_cpu_medium"]
    resources:
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*4:02}:59:59',
        bonus = 0
    shell:
        "bwa mem -t {threads} -R \"@RG\\tID:{wildcards.sample}_shortreads\\tSM:{wildcards.sample}\" "
            " {params.prefix} {input.reads} "
            "| samtools view -u -F 1792 "
            "| samtools sort -l 9 -m 2048M --threads {threads} -o {output.bam} --write-index /dev/stdin "
