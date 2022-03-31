
MINIMAP_CTG_REF_BASE = 'minimap2 -t {threads} -x asm20 -Y {params.sec_aln} '

MINIMAP_CTG_REF_BAM = MINIMAP_CTG_REF_BASE + '-a -L --MD --eqx '
MINIMAP_CTG_REF_BAM += '-R "@RG\\tID:{wildcards.sample}_{wildcards.hifi_type}_{wildcards.chrom}\\tSM:{wildcards.sample}" '
MINIMAP_CTG_REF_BAM += '{input.ref} {input.ctg} | samtools view -u -F 4 | '
MINIMAP_CTG_REF_BAM += 'samtools sort --threads {threads} -m {resources.sort_mem}M -l 9 -o {output} /dev/stdin '

MINIMAP_CTG_REF_PAF = MINIMAP_CTG_REF_BASE + '--cs -c --paf-no-hit {input.ref} {input.ctg} | pigz -p {threads} --best > {output}'

rule align_contigs_to_reference_paf:
    input:
        ctg = select_whole_genome_assembly,
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/{sub_folder}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/{sub_folder}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.mmap.rsrc'
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:59:00',
    params:
        sec_aln = "-p 0.95 --secondary=yes -N 1"
    shell:
        MINIMAP_CTG_REF_PAF


rule align_contigs_to_reference_bam:
    input:
        ctg = 'output/hybrid/renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam.mmap.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 32768 + 24576 * attempt,
        walltime = lambda wildcards, attempt: f'{23 * attempt}:59:00',
        sort_mem = lambda wildcards, attempt: 1024 * attempt
    params:
        sec_aln = "-p 0.95 --secondary=yes -N 1"
    shell:
        MINIMAP_CTG_REF_BAM


rule aggregate_contig_alignments:
    """
    This rule is essentially a preprocessing rule for subsetting
    the whole-genome assemblies to just chrY. This rule computes
    some basic statistics for all assembly contigs.
    """
    input:
        paf = 'output/alignments/contigs-to-ref/{sub_folder}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    output:
        tsv = 'output/alignments/contigs-to-ref/{sub_folder}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.ctg-agg.tsv'
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    conda:
        '../envs/pyscript.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = find_script_path('aggregate_contig_aln.py')
    shell:
        '{params.script_exec} --paf {input.paf} --output {output.tsv}'


rule align_reference_seq_classes:
    input:
        ref = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        ctg = 'references_derived/{seq_classes}.fasta'
    output:
        paf = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt}:59:00',
    params:
        sec_aln = "--secondary=yes -N 10000"
    shell:
        MINIMAP_CTG_REF_PAF
