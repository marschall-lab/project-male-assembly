
MINIMAP_CTG_REF_BASE = 'minimap2 -t {threads} -x asm20 -Y -m 10000 --end-bonus=100 -p 0.95 --secondary=yes -N 1 '

MINIMAP_CTG_REF_BAM = MINIMAP_CTG_REF_BASE + '-a -L --MD --eqx '
MINIMAP_CTG_REF_BAM += '-R "@RG\\tID:{wildcards.sample}_{wildcards.hifi_type}_{wildcards.chrom}\\tSM:{wildcards.sample}" '
MINIMAP_CTG_REF_BAM += '{input.ref} {input.ctg} | samtools view -u -F 4 | '
MINIMAP_CTG_REF_BAM += 'samtools sort --threads {threads} -l 9 -o {output} /dev/stdin '

MINIMAP_CTG_REF_PAF = MINIMAP_CTG_REF_BASE + '--cs -c --paf-no-hit {input.ref} {input.ctg} | pigz -p {threads} --best > {output}'

rule align_contigs_to_reference_paf:
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.mmap.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 20480 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:59:00',
    shell:
        MINIMAP_CTG_REF_PAF


rule align_contigs_to_reference_bam:
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.sort.bam'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam.mmap.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 20480 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:59:00',
    shell:
        MINIMAP_CTG_REF_BAM
