
rule call_variants_deepvariant_assm:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        fai = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.fai',
        bam = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam.bai',
    output:
        vcf = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.dv.vcf.gz',
        tbi = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.dv.vcf.gz.tbi',
        gvcf = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.dv.g.vcf.gz',
        tbi2 = 'output/variant_calls/00_dv_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.dv.g.vcf.gz.tbi',
    log:
        'log/output/variant_calls/00_dv_{other_reads}/{sample}.{hifi_type}.{ont_type}.na.chrY.deepvariant.log'
    benchmark:
        'rsrc/output/variant_calls/00_dv_{other_reads}/{sample}.{hifi_type}.{ont_type}.na.chrY.deepvariant.rsrc'
    wildcard_constraints:
        other_reads = 'HIFIRW'
    singularity:
        f"{config['container_store']}/{config['container']['deepvariant']}"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt:02}:59:59'
    params:
        model = 'PACBIO'
    shell:
        '/opt/deepvariant/bin/run_deepvariant --model_type="{params.model}" --ref={input.assm} '
            '--num_shards={threads} --reads={input.bam} --sample_name={wildcards.sample} '
            '--output_vcf={output.vcf} --output_gvcf={output.gvcf} '
            '--intermediate_results_dir=$TMPDIR &> {log} '


rule call_variants_pepper_assm:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        fai = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.fai',
        bam = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam.bai',
    output:
        vcf = 'output/variant_calls/00_pr_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.pr.vcf.gz',
        tbi = 'output/variant_calls/00_pr_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.pr.vcf.gz.tbi',
        gvcf = 'output/variant_calls/00_pr_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.pr.g.vcf.gz',
        tbi2 = 'output/variant_calls/00_pr_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.pr.g.vcf.gz.tbi',
    log:
        'log/output/variant_calls/00_pr_{other_reads}/{sample}.{hifi_type}.{ont_type}.na.chrY.pepper.log'
    benchmark:
        'rsrc/output/variant_calls/00_pr_{other_reads}/{sample}.{hifi_type}.{ont_type}.na.chrY.pepper.rsrc'
    wildcard_constraints:
        other_reads = 'ONTUL'
    singularity:
        f"{config['container_store']}/{config['container']['pepperdv']}"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*2:02}:59:59'
    params:
        outdir = lambda wildcards, output: pathlib.Path(output.vcf).parent,
        prefix = lambda wildcards, output: pathlib.Path(output.vcf).name.rsplit('.', 2)[0],
        preset = 'ont_r9_guppy5_sup'
    shell:
        'run_pepper_margin_deepvariant call_variant --bam {input.bam} --fasta {input.assm} --output_dir {params.outdir} '
            '--threads {threads} --{params.preset} --sample_name {wildcards.sample} --output_prefix {params.prefix} '
            '--skip_final_phased_bam --gvcf &> {log}'


rule select_hiq_het_snvs:
    input:
        vcf = 'output/variant_calls/00_{caller}_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.{caller}.vcf.gz',
        tbi = 'output/variant_calls/00_{caller}_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.{caller}.vcf.gz.tbi'
    output:
        vcf = 'output/variant_calls/10_filter_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.{caller}-HET-SNV.vcf.gz',
        tbi = 'output/variant_calls/10_filter_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.{caller}-HET-SNV.vcf.gz.tbi',
        stats = 'output/variant_calls/10_filter_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.{caller}-HET-SNV.stats',
    wildcard_constraints:
        caller = '(pr|dv)'
    conda:
        '../envs/biotools.yaml'
    params:
        min_qual = 10
    shell:
        'bcftools filter --include "QUAL>={params.min_qual}" -Ov {input.vcf} | bcftools view -Oz --types snps --genotype het > {output.vcf}'
            ' && '
        'tabix --preset vcf {output.vcf}'
            ' && '
        'bcftools stats {output.vcf} > {output.stats}'


rule convert_het_snv_to_tsv:
    input:
        vcf = 'output/variant_calls/10_filter_{other_reads}/{sample}/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{caller}-HET-SNV.vcf.gz',
    output:
        tsv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.{caller}-errors.tsv'
    run:
        import gzip
        import pandas as pd
        source = 'DeepVariant' if wildcards.caller == 'dv' else 'PEPPER'

        rows = []
        with gzip.open(input.vcf, 'rt') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                row = {
                    'chrom': parts[0],
                    'end': int(parts[1]),
                    'start': int(parts[1]) - 1,
                    'est_size': 1,
                    'sample': wildcards.sample,
                    'source': source,
                    'reads': wildcards.other_reads
                }
                rows.append(row)
        col_names = ['chrom', 'start', 'end', 'est_size', 'sample', 'source', 'reads']
        df = pd.DataFrame.from_records(rows, columns=col_names)
        df.sort_values(['chrom', 'start'], ascending=True, inplace=True)
        df.to_csv(output.tsv, sep='\t', header=True, index=False)
    # END OF RUN BLOCK

