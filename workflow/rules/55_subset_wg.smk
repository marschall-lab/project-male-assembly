

rule determine_chry_contigs:
    """
    This makes only sense relative to a complete chrY, i.e. T2T,
    hence reference is hard-coded
    """
    input:
        agg_ctg_aln = expand(
            'output/alignments/contigs-to-ref/{{sample_info}}_{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}_aln-to_{reference}.ctg-agg.tsv',
            reference='T2TXY',
            chrom='wg',
            mapq='na'
        ),
        agg_motifs = expand(
            'output/motif_search/20_target_agg/{{sample_info}}_{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}.{motif}.agg-trg.tsv',
            motif=config['motif_search'],
            chrom='wg',
            mapq='na'
        )
    output:
        table = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.stats.tsv',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        bed = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.bed',
    conda:
        '../envs/pyscript.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = find_script_path('determine_y_contigs.py')
    shell:
        '{params.script_exec} --agg-align {input.agg_ctg_aln} --agg-motif {input.agg_motifs} '
            '--out-stats {output.table} --out-names {output.names} --out-bed {output.bed}'


rule extract_y_contigs:
    input:
        wg_assm = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg/assembly.fasta',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
    output:
        sub_assm = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'seqtk subseq {input.wg_assm} {input.names} > {output.sub_assm}'


rule extract_contig_alignments_paf:
    input:
        paf = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.paf.gz',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
    output:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -F -f {input.names} {input.paf} | gzip > {output.paf}'


rule extract_contig_alignments_bam:
    """
    TODO: check if sorted input generates sorted output
    """
    input:
        bam = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.sort.bam',
        bai = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.sort.bam.bai',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
    output:
        bam = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.sort.bam'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt
    shell:
        'samtools view -u --qname-file {input.names} {input.bam} | samtools sort -m 2048M -l 9 -@ 2 -o {output.bam} /dev/stdin'


rule extract_read_alignments_paf:
    input:
        paf = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.paf.gz',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
    output:
        paf = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -F -f {input.names} {input.paf} | gzip > {output.paf}'


rule extract_read_alignments_bam:
    """
    TODO: check if sorted input generates sorted output
    """
    input:
        bam = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai',
        bed = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.bed',
    output:
        bam = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt ** attempt:02}:59:00',
        sort_mem = lambda wildcards, attempt: 4096 * attempt
    shell:
        'samtools view -u -L {input.bed} {input.bam} | samtools sort -m {resources.sort_mem}M -l 9 -@ {threads} -o {output.bam} /dev/stdin'


rule subset_motif_hits:
    input:
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        bed = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm-hiq.bed',
    output:
        bed = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
    shell:
        'grep -F -f {input.names} {input.bed} | sort -V -k1 -k2n,3n > {output.bed}'


rule extract_motif_hit_sequences:
    input:
        fasta = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        bed = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
    output:
        fasta = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.hiq-seq.fasta',
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq {input.fasta} {input.bed} > {output.fasta}'
