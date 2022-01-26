
MINIMAP_CTG_REF_BASE = 'minimap2 -t {threads} -x asm20 -Y -m 10000 --end-bonus=100 -p 0.95 --secondary=yes -N 1 '

MINIMAP_CTG_REF_BAM = MINIMAP_CTG_REF_BASE + '-a -L --MD --cs '
MINIMAP_CTG_REF_BAM += '-R "@RG\\tID:{wildcards.sample}_{wildcards.hifi_type}_{wildcards.chrom}\\tSM:{wildcards.sample}" '
MINIMAP_CTG_REF_BAM += '{input.ref} {input.ctg} | samtools sort --threads {threads} -l 9 -o {output} /dev/stdin '

MINIMAP_CTG_REF_PAF = MINIMAP_CTG_REF_BASE + '-c --paf-no-hit {input.ref} {input.ctg} | pigz -p {threads} --best > {output}'

MINIMAP_CTG_REF_CALLS = {
    'sort.bam': MINIMAP_CTG_REF_BAM,
    'paf.gz': MINIMAP_CTG_REF_PAF
}

rule align_contigs_to_reference:
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.{ext}'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.{ext}.mmap.rsrc'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:59:00',
    params:
        mm_call = lambda wildcards: MINIMAP_CTG_REF_CALLS[wildcards.ext]
    shell:
        '{params.mm_call}'
