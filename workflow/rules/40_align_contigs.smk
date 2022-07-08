
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
        aln = 'output/alignments/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.mmap.rsrc'
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


rule align_contigs_close_samples:
    input:
        na19317 = 'output/subset_wg/20_extract_contigs/NA19317.{hifi_type}.{ont_type}.na.chrY.fasta',
        na19347 = 'output/subset_wg/20_extract_contigs/NA19347.{hifi_type}.{ont_type}.na.chrY.fasta',
    output:
        ref_317 = 'output/alignments/contigs-to-contigs/NA19347.{hifi_type}.{ont_type}.na.chrY_aln-to_NA19317.paf.gz',
        ref_347 = 'output/alignments/contigs-to-contigs/NA19317.{hifi_type}.{ont_type}.na.chrY_aln-to_NA19347.paf.gz',
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:59:00',
    params:
        sec_aln = "-p 0.95 --secondary=yes -N 1"
    shell:
        'minimap2 -t {threads} -x asm20 -Y {params.sec_aln} '
        '--cs -c --paf-no-hit {input.na19317} {input.na19347} | pigz -p {threads} --best > {output.ref_317}'
            ' && '
        'minimap2 -t {threads} -x asm20 -Y {params.sec_aln} '
        '--cs -c --paf-no-hit {input.na19347} {input.na19317} | pigz -p {threads} --best > {output.ref_347}'


rule cache_close_contig_alignments:
    """
    Store all contig alignments
    for later plotting (more efficient binning);
    for the two benchmark samples, keep primary
    and secondary alignment information intact
    """
    input:
        paf = 'output/alignments/contigs-to-contigs/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{query_sample}.paf.gz'
    output:
        hdf = 'output/alignments/contigs-to-contigs/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{query_sample}.cache.h5'
    run:
        import pandas as pd
        import numpy as np

        PAF_COLUMN_NAMES = [
            'query_chrom', 'query_length', 'query_start', 'query_end', 'orientation',
            'target_chrom', 'target_length', 'target_start', 'target_end', 'res_matches', 'block_length',
            'mapq', 'tag_nm', 'tag_ms', 'tag_as', 'tag_nn', 'tag_tp'
        ]
        PAF_USE_COLS = list(range(0,17))
        assert len(PAF_COLUMN_NAMES) == len(PAF_USE_COLS)

        df = pd.read_csv(input.paf, sep='\t', header=None, names=PAF_COLUMN_NAMES, usecols=PAF_USE_COLS)
        df['tag_tp'] = df['tag_tp'].str.lower()
        df.sort_values(['target_name', 'target_start'], ascending=True, inplace=True)

        with pd.HDFStore(output.hdf, mode='w', complib='blosc', complevel=9) as hdf:
            pass

        target_contigs = []
        for ctg_order_key, (target_ctg, ctg_aln) in enumerate(df.groupby('target_chrom'), start=1):
            ctg_key = f'{ctg_order_key:03}'
            ctg_size = ctg_aln.at[ctg_aln.index[0], 'target_length']

            target_contigs.append((ctg_key, target_ctg, ctg_size))

            for aln_type, type_label in zip([['tp:a:p', 'tp:a:i'], ['tp:a:s']], ['primary', 'secondary']):
                sub_aln = ctg_aln.loc[ctg_aln['tag_tp'].isin(aln_type), :]
                aln_cov = np.zeros(ctg_size, dtype=np.int8)
                for tstart, tend in sub_aln[['target_start', 'target_end']].itertuples(index=False):
                    aln_cov[tstart:tend] += 1
                    if aln_cov.max() > 254:
                        raise
                
                with pd.HDFStore(output.hdf, mode='a', complib='blosc', complevel=9) as hdf:
                    store_key = f'{wildcards.sample}/{ctg_key}/{type_label}'
                    hdf.put(store_key, pd.Series(aln_cov), format='fixed')

        target_contigs = pd.DataFrame.from_records(target_contigs, names=['order_key', 'contig_name', 'contig_size'])
        with pd.HDFStore(output.hdf, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put('contigs', target_contigs, format='fixed')

    # END OF RUN BLOCK


rule align_contigs_to_reference_bam:
    input:
        ctg = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
        ref = lambda wildcards: select_reference_genome(wildcards.reference)
    output:
        aln = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam'
    benchmark:
        'rsrc/output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam.mmap.rsrc'
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
        paf = 'output/alignments/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    output:
        tsv = 'output/alignments/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.ctg-agg.tsv'
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
        ref = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        ctg = 'references_derived/{seq_classes}.fasta'
    output:
        paf = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.paf.gz'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 + 4096 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
    params:
        sec_aln = "--secondary=yes -N 10000"
    shell:
        MINIMAP_CTG_REF_PAF


rule convert_seq_class_alignments_to_bed:
    input:
        'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.paf.gz'
    output:
        'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.{alignments}.bed'
    wildcard_constraints:
        alignments = '(primary|secondary)'
    params:
        select_aln = lambda wildcards: 'tp:A:P' if wildcards.alignments == 'primary' else 'tp:A:S'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
    shell:
        'zgrep -F "{params.select_aln}" {input} | cut -f 1-12 | '
        'awk \'BEGIN{{OFS = "\\t"}} {{print $6,$8,$9,$1,$12,$5,$3,$4,$2,"{wildcards.alignments}"}}\' | '
        ' sort -k1 -k2,3n > {output}'


rule merge_seq_class_bed_files:
    input:
        pri = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.primary.bed',
        sec = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.secondary.bed'
    output:
        mrg = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.na.chrY.merged.bed'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
    shell:
        'cat {input.pri} {input.sec} | sort -k1 -k2,3n > {output.mrg}'
