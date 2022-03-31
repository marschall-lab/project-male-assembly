

rule determine_chry_contigs:
    """
    This makes only sense relative to a complete chrY, i.e. T2T,
    hence reference is hard-coded

    NB: this rule is only executed on the raw Verkko output, i.e.,
    it it a prerequisite for identifying and renaming the assembled
    chrY contigs, but its output is not used afterwards
    """
    input:
        agg_ctg_aln = expand(
            'output/alignments/contigs-to-ref/00_raw/{{sample_info}}_{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}_aln-to_{reference}.ctg-agg.tsv',
            reference='T2TXY',
            chrom='wg',
            mapq='na'
        ),
        agg_motifs = expand(
            'output/motif_search/20_target_agg/00_raw/{{sample_info}}_{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}.{motif}.agg-trg.tsv',
            motif=config['contig_id_motifs'],
            chrom='wg',
            mapq='na'
        )
    output:
        table = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.stats.tsv',
        names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        bed = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
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
        ctg_names = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        ctg_aln = 'output/alignments/contigs-to-ref/00_raw/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.paf.gz',
        seq_classes = lambda wildcards: f'references_derived/{config["reference_y_seq_classes"]["T2TXY"]}.bed',
        ref_cov = 'output/eval/contigs-to-ref/00_raw/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.ref-cov.tsv',
    output:
        new_names = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        name_maps = multiext(
            'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names',
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
        wg_hifi_cov = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg/assembly.hifi-coverage.csv',
        wg_ont_cov = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg/assembly.ont-coverage.csv',
        ren_json = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.json',
        ren_sed = 'output/subset_wg/15_order_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.sed',
        sub_bed = 'output/subset_wg/10_find_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
    output:
        ren_assm = 'output/hybrid/renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
        ren_hifi_cov = 'output/hybrid/renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.hifi-coverage.csv',
        ren_ont_cov = 'output/hybrid/renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.ont-coverage.csv',
        sub_assm = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        sub_bed = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
        sub_names = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    conda:
        '../envs/pyscript.yaml'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt}:59:59',
        mem_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        script_exec = find_script_path('rename_extract_assembly.py')
    shell:
        '{params.script_exec} --input-fasta {input.wg_assm} --name-map {input.ren_json} '
        '--out-wg {output.ren_assm} --out-sub {output.sub_assm}'
            ' && '
        'sed -f {input.ren_sed} {input.wg_ont_cov} > {output.ren_ont_cov}'
            ' && '
        'sed -f {input.ren_sed} {input.wg_hifi_cov} > {output.ren_hifi_cov}'
            ' && '
        'sed -f {input.ren_sed} {input.sub_bed} > {output.sub_bed}'
            ' && '
        'cut -f 1 {output.sub_bed} | sort > {output.sub_names}'


rule extract_contig_alignments_paf:
    input:
        paf = 'output/alignments/contigs-to-ref/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.paf.gz',
        names = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    output:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | sort -k1 -k3n,4n | gzip > {output.paf}'


rule extract_contig_alignments_bam:
    """
    TODO: check if sorted input generates sorted output

    Renaming: de novo chrY contigs are the queries, i.e.
    need to dump to text/SAM, rename, and then recompress
    to BAM.
    """
    input:
        bam = 'output/alignments/contigs-to-ref/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.bam',
        bai = 'output/alignments/contigs-to-ref/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.bam.bai',
        names = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    output:
        bam = 'output/subset_wg/30_extract_ctgaln/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.bam'
    conda:
        '../envs/biotools.yaml'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        walltime = lambda wildcards, attempt: f'{6 * attempt}:59:00'
    shell:
        'samtools view -u --qname-file {input.names} {input.bam} | '
        'samtools sort -m 2048M -l 9 --output-fmt BAM --threads 4 -o {output.bam} /dev/stdin'


rule extract_read_alignments_paf:
    input:
        paf = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.paf.gz',
        names = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    output:
        paf = 'output/subset_wg/40_extract_rdaln/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | pigz -p 4 --best > {output.paf}'


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
        bed = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
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
        'samtools view -u -L {input.bed} {input.bam} | samtools sort -m {resources.sort_mem}M -l 4 -@ {threads} -o {output.bam} /dev/stdin'


rule subset_motif_hits:
    """
    If the input BED is empty (no hits above threshold),
    touch the output file.
    The TSV also contains low-scoring hits, so it is unlikely
    to be empty, but better safe than sorry.
    """
    input:
        names = 'output/subset_wg/20_extract_contigs/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        tsv = 'output/motif_search/10_norm/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm.tsv',
        bed = 'output/motif_search/10_norm/10_renamed/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm-hiq.bed',
    output:
        tsv = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm.tsv',
        bed = 'output/subset_wg/50_subset_motif/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
    shell:
        'if [ -s {input.bed} ]; then '
        'grep -w -F -f {input.names} {input.bed} | sort -V -k1 -k2n,3n > {output.bed} ; '
        'else '
        'touch {output.bed} ; '
        'fi ;'
        'if [ -s {input.tsv} ]; then '
        'grep -w -F -f {input.names} {input.tsv} | sort -V -k1 -k5n,6n > {output.tsv} ; '
        'else '
        'touch {output.tsv} ; '
        'fi ;'


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
        'zgrep -E "{params.chroms}" {input.paf} | pigz -p 2 --best > {output.paf}'
