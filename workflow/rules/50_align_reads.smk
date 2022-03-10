
MINIMAP_READ_ASSM_BASE = 'minimap2 -t {threads} -x {params.preset} -Y -p 0.95 --secondary=yes -N 1 --cap-kalloc=1g -K4g -I8g '

MINIMAP_READ_ASSM_BAM = MINIMAP_READ_ASSM_BASE + '-a -L --MD --eqx '
MINIMAP_READ_ASSM_BAM += '-R "@RG\\tID:{wildcards.sample}_{wildcards.hifi_type}_{wildcards.chrom}_{wildcards.other_reads}\\tSM:{wildcards.sample}" '
MINIMAP_READ_ASSM_BAM += '{input.ctg} {input.reads} | samtools view -u -F 4 | '
MINIMAP_READ_ASSM_BAM += 'samtools sort --threads {threads} -m {resources.sort_mem}M -l 9 -o {output} /dev/stdin'

MINIMAP_READ_ASSM_PAF = MINIMAP_READ_ASSM_BASE + '--cs -c --paf-no-hit {input.ctg} {input.reads} | pigz -p {threads} --best > {output}'

MINIMAP_PRESETS = {
    'HIFIRW': 'map-hifi',
    'HIFIAF': 'map-hifi',
    'HIFIEC': 'map-hifi',
    'ONTUL': 'map-ont',
    'ONTEC': 'map-hifi',
    'OHEC': 'map-hifi'
}


def select_input_reads(wildcards):

    if wildcards.other_reads == 'HIFIEC':
        reads = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.{chrom}/hifi-corrected.fasta'.format(**dict(wildcards))
    else:
        reads = SAMPLE_INFOS[wildcards.sample][wildcards.other_reads]
    return reads


rule align_reads_to_assembly_paf:
    """
    This rule only aligns whole-genome read sets
    """
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.{chrom}/assembly.fasta',
        reads = select_input_reads,
    output:
        'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.paf.gz'
    log:
        'log/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.mmap-paf.log'
    benchmark:
        'rsrc/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.mmap-paf.rsrc'
    wildcard_constraints:
        chrom = 'wg'
    threads: config['num_cpu_medium']
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 65536 + 48576 * attempt,
        walltime = lambda wildcards, attempt: f'{47*attempt}:59:00'
    params:
        preset = lambda wildcards: MINIMAP_PRESETS[wildcards.other_reads]
    shell:
        MINIMAP_READ_ASSM_PAF


rule align_reads_to_assembly_bam:
    """
    This rule only aligns whole-genome read sets
    """
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.{chrom}/assembly.fasta',
        reads = select_input_reads,
    output:
        'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.bam'
    log:
        'log/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.mmap-bam.log'
    benchmark:
        'rsrc/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.mmap-bam.rsrc'
    wildcard_constraints:
        chrom = 'wg'
    threads: config['num_cpu_medium']
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 65536 + 48576 * attempt,
        walltime = lambda wildcards, attempt: f'{47*attempt}:59:00',
        sort_mem = lambda wildcards, attempt: 2048 * attempt
    params:
        preset = lambda wildcards: MINIMAP_PRESETS[wildcards.other_reads]
    shell:
        MINIMAP_READ_ASSM_BAM


rule align_subset_reads_to_assembly_paf:
    """
    This rule only aligns read subsets
    to a chromosome assembly (Y, potentially, X and XY),
    but for all possible combinations. The most important
    one is: the assembly reads, i.e. HIFIAF (by default)
    and ONTUL
    """
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        reads = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{other_reads}.{chrom}-reads.{mapq}.fasta.gz',
    output:
        'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.paf.gz'
    log:
        'log/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.mmap-paf.log'
    benchmark:
        'rsrc/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.mmap-paf.rsrc'
    wildcard_constraints:
        mapq = 'mq00',
        chrom = '(chrY|chrX|chrXY)'
    threads: config['num_cpu_low']
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    params:
        preset = lambda wildcards: MINIMAP_PRESETS[wildcards.other_reads]
    shell:
        MINIMAP_READ_ASSM_PAF


rule align_subset_reads_to_assembly_bam:
    """
    This rule only aligns read subsets
    to a chromosome assembly (Y, potentially, X and XY),
    but for all possible combinations. The most important
    one is: the assembly reads, i.e. HIFIAF (by default)
    and ONTUL
    """
    input:
        ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        reads = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{other_reads}.{chrom}-reads.{mapq}.fasta.gz',
    output:
        'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.bam'
    log:
        'log/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.mmap-bam.log'
    benchmark:
        'rsrc/output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.mmap-bam.rsrc'
    wildcard_constraints:
        mapq = 'mq00',
        chrom = '(chrY|chrX|chrXY)'
    threads: config['num_cpu_low']
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    params:
        preset = lambda wildcards: MINIMAP_PRESETS[wildcards.other_reads]
    shell:
        MINIMAP_READ_ASSM_BAM
