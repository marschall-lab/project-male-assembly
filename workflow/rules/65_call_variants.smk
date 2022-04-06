
rule call_variants_deepvariant_assm:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        fai = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.fai',
        bam = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam.bai',
    output:
        vcf = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.dv.vcf.gz',
        gvcf = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.dv.gvcf.gz',
    log:
        'log/output/variant_calls/00_dv_{other_reads}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.deepvariant.log'
    benchmark:
        'rsrc/output/variant_calls/00_dv_{other_reads}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.deepvariant.rsrc'
    wildcard_constraints:
        other_reads = 'HIFIRW'
    singularity:
        f'{config['container_store']}/{config['deepvariant']}'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 24576 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt**3}:59:59'
    params:
        model = 'PACBIO'
    shell:
        '/opt/deepvariant/bin/run_deepvariant --model_type="{params.model}" --ref={input.assm} '
            '--num_shards={threads} --reads={input.bam} --sample_name={wildcards.sample} '
            '--output_vcf={output.vcf} --output_gvcf={output.gvcf} '
            '--intermediate_results_dir=$TMPDIR &> {log} '


rule call_variants_pepper_assm:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        fai = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.fai',
        bam = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam.bai',
    output:
        outdir = directory('output/variant_calls/00_pr_{other_reads}/{sample}')
    log:
        'log/output/variant_calls/00_pr_{other_reads}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.pepper.log'
    benchmark:
        'rsrc/output/variant_calls/00_pr_{other_reads}/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.pepper.rsrc'
    wildcard_constraints:
        other_reads = 'ONTUL'
    singularity:
        f'{config['container_store']}/{config['pepperdv']}'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 24576 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt**3}:59:59'
    params:
        prefix = lambda wildcards: f'{wildcards.sample_info}_{wildcards.sample}.{wildcards.hifi_type}.{wildcards.ont_type}.na.chrY.pr',
        preset = 'ont_r9_guppy5_sup'
    shell:
        'run_pepper_margin_deepvariant call_variant --bam {input.bam} --fasta {input.assm} --output_dir {output.outdir} '
            '--threads {threads} --{params.preset} --sample_name {wildcards.sample} --output_prefix {params.prefix} '
            '--skip_final_phased_bam --gvcf &> {log}'
