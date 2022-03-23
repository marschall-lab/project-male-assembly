

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
            motif=config['contig_id_motifs'],
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


rule determine_contig_order:
    """
    Same as above: ordering contigs relative to a reference assembly
    makes only sense for a complete assembly, i.e. T2T
    """
    input:
        ctg_names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        ctg_aln = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.paf.gz',
        seq_classes = lambda wildcards: f'references_derived/{config["reference_y_seq_classes"]["T2TXY"]}.bed',
        ref_cov = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.ref-cov.tsv',
    output:
        new_names = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        name_maps = multiext(
            'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names',
            '.otn-map.tsv', '.otn-map.json', '.otn-map.sed',
            '.nto-map.tsv', '.nto-map.json', '.nto-map.sed'
        )
    conda:
        '../envs/pyscript.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = find_script_path('determine_contig_order.py')
    shell:
        '{params.script_exec} --sample-name {wildcards.sample} --names {input.ctg_names} '
            '--seq-classes {input.seq_classes} --class-coverage {input.ref_cov} '
            '--paf {input.ctg_aln} --output {output.new_names} '
            '--dump-mappings tsv json sed'


rule extract_y_contigs:
    input:
        wg_assm = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg/assembly.fasta',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        sub_assm = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'seqtk subseq {input.wg_assm} {input.names} | sed -f {input.rename} > {output.sub_assm}'


rule extract_contig_alignments_paf:
    input:
        paf = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.paf.gz',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | sed -f {input.rename} | sort -k1 -k3n,4n | gzip > {output.paf}'


rule extract_contig_alignments_bam:
    """
    TODO: check if sorted input generates sorted output

    Renaming: de novo chrY contigs are the queries, i.e.
    need to dump to text/SAM, rename, and then recompress
    to BAM.
    """
    input:
        bam = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.sort.bam',
        bai = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.sort.bam.bai',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        bam = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.sort.bam'
    conda:
        '../envs/biotools.yaml'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt
    shell:
        'samtools view --with-header --qname-file {input.names} {input.bam} | '
        'sed -f {input.rename} | '
        'samtools sort -m 2048M -l 9 --output-fmt BAM --threads 4 -o {output.bam} /dev/stdin'


rule extract_read_alignments_paf:
    input:
        paf = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.paf.gz',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        paf = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | sed -f {input.rename} | pigz -p 4 --best > {output.paf}'


rule extract_read_alignments_bam:
    """
    TODO: check if sorted input generates sorted output

    Renaming: de novo chrY contigs are the targets, i.e.
    need to reheader the output BAM file after extracting
    the read alignments.
    """
    input:
        bam = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai',
        bed = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.bed',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        tmp_bam = temp('output/subset_wg/40_extract_rdaln/tmp/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam'),
        bam = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt ** attempt:02}:59:00',
        sort_mem = lambda wildcards, attempt: 4096 * attempt
    shell:
        'samtools view -u -L {input.bed} {input.bam} | samtools sort -m {resources.sort_mem}M -l 4 -@ {threads} -o {output.tmp_bam} /dev/stdin'
            ' && '
        'samtools view -H {output.tmp_bam} | sed -f {input.rename} | samtools reheader /dev/stdin {output.tmp_bam} > {output.bam}'


rule subset_motif_hits:
    """
    If the input BED is empty (no hits above threshold),
    touch the output file
    """
    input:
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.txt',
        bed = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm-hiq.bed',
        rename = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.chrY.names.otn-map.sed',
    output:
        bed = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
    shell:
        'if [ -s {input.bed} ]; then '
        'grep -w -F -f {input.names} {input.bed} | sed -f {input.rename} | sort -V -k1 -k2n,3n > {output.bed} ; '
        'else '
        'touch {output.bed} ; '
        'fi'


rule extract_motif_hit_sequences:
    """
    NB: all input files to this rule have been renamed before, i.e.,
    no need to rename chrY contigs here

    NB: if the input BED is empty, seqtk will simply produce an empty output FASTA
    """
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


REF_CHRY = {
    'GRCh38': ['chrY', 'chrY_KI270740v1_random'],
    'T2TXY': ['chrY']
}

rule extract_read_ref_alignments_bam:
    input:
        bam = 'output/alignments/reads-to-ref/{sample_info}_{sample}.{other_reads}_aln-to_{reference}.bam',
        bai = 'output/alignments/reads-to-ref/{sample_info}_{sample}.{other_reads}_aln-to_{reference}.bam.bai',
    output:
        bam = 'output/subset_wg/60_subset_rdref/{sample_info}_{sample}.{other_reads}_aln-to_{reference}.chrY.bam'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        chroms = lambda wildcards: ' '.join(REF_CHRY[wildcards.reference])
    shell:
        'samtools view -b {input.bam} {params.chroms} > {output.bam}'


rule extract_read_ref_alignments_paf:
    input:
        paf = 'output/alignments/reads-to-ref/{sample_info}_{sample}.{other_reads}_aln-to_{reference}.paf.gz',
    output:
        paf = 'output/subset_wg/60_subset_rdref/{sample_info}_{sample}.{other_reads}_aln-to_{reference}.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        chroms = lambda wildcards: '(' + '|'.join(REF_CHRY[wildcards.reference]) + ')'
    shell:
        'zgrep -E {params.chroms} {input.paf} | pigz -p 2 --best > {output.paf}'
