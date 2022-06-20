

rule determine_chrom_contigs:
    """
    This makes only sense relative to a complete chrY, i.e. T2T,
    hence reference is hard-coded

    NB: this rule is only executed on the raw Verkko output, i.e.,
    it it a prerequisite for identifying and renaming the assembled
    chrY contigs, but its output is not used afterwards

    Update: 2022-05-27
    chrX also needs to be extracted to take a look at the PAR regions.
    However, since no sequence motifs are available to aid in the
    identification of the respective contigs, the "agg_motifs"
    input will effectively be ignored by the script, and no
    further interpretation of the contigs is performed as part
    of this pipeline (all downstream).

    """
    input:
        agg_ctg_aln = expand(
            'output/alignments/contigs-to-ref/00_raw/{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}_aln-to_{reference}.ctg-agg.tsv',
            reference='T2TXY',
            chrom='wg',
            mapq='na'
        ),
        agg_motifs = expand(
            'output/motif_search/20_target_agg/00_raw/{{sample}}.{{hifi_type}}.{{ont_type}}.{mapq}.{chrom}.{motif}.agg-trg.tsv',
            motif=config['contig_id_motifs'],
            chrom='wg',
            mapq='na'
        )
    output:
        table = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.stats.tsv',
        names = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt',
        bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
    wildcard_constraints:
        chrom = '(chrY|chrX)'
    conda:
        '../envs/pyscript.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = find_script_path('identify_contigs.py')
    shell:
        '{params.script_exec} --select-chrom {wildcards.chrom} --agg-align {input.agg_ctg_aln} '
            '--agg-motif {input.agg_motifs} --out-stats {output.table} --out-names {output.names} --out-bed {output.bed}'


rule determine_contig_order:
    """
    Same as above: ordering contigs relative to a reference assembly
    makes only sense for a complete assembly, i.e. T2T
    """
    input:
        ctg_names = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt',
        ctg_aln = 'output/alignments/contigs-to-ref/00_raw/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.paf.gz',
        seq_classes = lambda wildcards: f'references_derived/{config["reference_y_seq_classes"]["T2TXY"]}.bed',
        ref_cov = 'output/eval/contigs-to-ref/00_raw/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_T2TXY.ref-cov.tsv',
    output:
        new_names = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt',
        name_maps = multiext(
            'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names',
            '.otn-map.tsv', '.otn-map.json', '.otn-map.sed',
            '.nto-map.tsv', '.nto-map.json', '.nto-map.sed'
        )
    wildcard_constraints:
        chrom = '(chrY|chrX)'
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
            '--dump-mappings tsv json sed --process-chrom {wildcards.chrom}'


rule extract_chrom_contigs:
    input:
        wg_assm = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.fasta',
        ren_y_json = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.json',
        ren_x_json = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.names.otn-map.json',
    output:
        ren_assm = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
        sub_y_assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        sub_x_assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.fasta',
    conda:
        '../envs/pyscript.yaml'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt}:59:59',
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        script_exec = find_script_path('rename_extract_assembly.py')
    shell:
        '{params.script_exec} --input-fasta {input.wg_assm} --out-wg {output.ren_assm} '
            '--name-map {input.ren_y_json} {input.ren_x_json} '
            '--out-sub {output.sub_y_assm} {output.sub_x_assm}'


rule rename_verkko_coverage_tables:
    input:
        wg_hifi_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.hifi-coverage.csv',
        wg_ont_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.ont-coverage.csv',
        ren_y_sed = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.sed',
        sub_y_bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
        ren_x_sed = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.names.otn-map.sed',
        sub_x_bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.bed',
    output:
        ren_hifi_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.hifi-coverage.csv',
        ren_ont_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.ont-coverage.csv',
        merge_sed = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY-chrX.names.otn-map.sed',
        sub_y_bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
        sub_y_names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        sub_x_bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.bed',
        sub_x_names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.names.txt'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt:02}:59:59',
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'cat {input.ren_y_sed} {input.ren_x_sed} > {output.merge_sed}'
            ' && '
        'sed -f {output.merge_sed} {input.wg_ont_cov} > {output.ren_ont_cov}'
            ' && '
        'sed -f {output.merge_sed} {input.wg_hifi_cov} > {output.ren_hifi_cov}'
            ' && '
        'sed -f {input.ren_y_sed} {input.sub_y_bed} > {output.sub_y_bed}'
            ' && '
        'cut -f 1 {output.sub_y_bed} | sort > {output.sub_y_names}'
            ' && '
        'sed -f {input.ren_x_sed} {input.sub_x_bed} > {output.sub_x_bed}'
            ' && '
        'cut -f 1 {output.sub_x_bed} | sort > {output.sub_x_names}'


rule extract_contig_alignments_paf:
    input:
        paf = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.paf.gz',
        names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    output:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.paf.gz'
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
        bam = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.bam',
        bai = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.bam.bai',
        names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt'
    output:
        bam = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{reference}.bam'
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
        paf = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.paf.gz',
        names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt'
    output:
        paf = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.paf.gz'
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    conda:
        '../envs/biotools.yaml'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | pigz -p 4 --best > {output.paf}'


rule extract_read_names:
    """
    This rule exists to run VerityMap, which does not support
    diploid genome assemblies, i.e., need to run only on
    chrY and chrX
    """
    input:
        paf = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.paf.gz'
    output:
        txt = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.read-names.txt'
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zcat {input.paf} | cut -f 1 | sort | uniq > {output.txt}'


rule extract_read_subset:
    input:
        reads = select_input_reads,
        names = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.read-names.txt'
    output:
        fasta = 'output/subset_wg/45_extract_reads/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.reads.fasta.gz'
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    conda:
        '../envs/biotools.yaml'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt
    shell:
        'zcat {input.reads} | seqtk subseq /dev/stdin {input.names} | seqtk seq -A -C | pigz --best -p {threads} > {output.fasta}'


rule extract_read_alignments_bam:
    """
    TODO: check if sorted input generates sorted output

    Renaming: de novo chrY contigs are the targets, i.e.
    need to reheader the output BAM file after extracting
    the read alignments.
    """
    input:
        bam = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai',
        bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
    output:
        bam = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
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
        names = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.names.txt',
        tsv = 'output/motif_search/10_norm/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm.tsv',
        bed = 'output/motif_search/10_norm/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg.{motif}.norm-hiq.bed',
    output:
        tsv = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm.tsv',
        bed = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
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
        fasta = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        bed = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.norm-hiq.bed',
    output:
        fasta = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.na.chrY.{motif}.hiq-seq.fasta',
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
        bam = 'output/alignments/reads-to-ref/{sample}.{other_reads}_aln-to_{reference}.bam',
        bai = 'output/alignments/reads-to-ref/{sample}.{other_reads}_aln-to_{reference}.bam.bai',
    output:
        bam = 'output/subset_wg/60_subset_rdref/{sample}.{other_reads}_aln-to_{reference}.chrY.bam'
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
        paf = 'output/alignments/reads-to-ref/{sample}.{other_reads}_aln-to_{reference}.paf.gz',
    output:
        paf = 'output/subset_wg/60_subset_rdref/{sample}.{other_reads}_aln-to_{reference}.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        chroms = lambda wildcards: '(' + '|'.join(REF_CHRY[wildcards.reference]) + ')'
    shell:
        'zgrep -E "{params.chroms}" {input.paf} | pigz -p 2 --best > {output.paf}'
