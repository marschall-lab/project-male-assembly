

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
            'output/alignments/contigs-to-ref/00_raw/{{sample}}.{{hifi_type}}.{{ont_type}}.{{mapq}}.{chrom}_aln-to_{reference}.ctg-agg.tsv',
            reference='T2TXY',
            chrom='wg',
        ),
        agg_motifs = expand(
            'output/motif_search/20_target_agg/00_raw/{{sample}}.{{hifi_type}}.{{ont_type}}.{{mapq}}.{chrom}.{motif}.agg-trg.tsv',
            motif=config['contig_id_motifs'],
            chrom='wg',
        )
    output:
        table = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.stats.tsv',
        names = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.txt',
        bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.bed',
    wildcard_constraints:
        chrom = '(chrY|chrX)',
        mapq = '(na|ha)'
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

    Update 2022-06-20
    - Edge case problem for samples NA19239 and HG03492 with manually renamed contigs:
    the renaming renders the SED/JSON/TSV name mapping generated by this rule invlid for one entry,
    which is only relevant for RepeatMasker (need temp copy of FASTA with original identifiers).
    - However, the output of this rule is needed for "extract_chrom_contigs", and can thus not be
    omitted.
    - Somewhat ugly but explicit solution: create the name mapping files after having
    dumped the subset FASTA files, and keep only the JSON output mapping of this rule
    """
    input:
        ctg_names = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.txt',
        ctg_aln = 'output/alignments/contigs-to-ref/00_raw/{sample}.{hifi_type}.{ont_type}.{mapq}.wg_aln-to_T2TXY.paf.gz',
        seq_classes = lambda wildcards: f'references_derived/{config["reference_y_seq_classes"]["T2TXY"]}.bed',
        ref_cov = 'output/eval/contigs-to-ref/00_raw/{sample}.{hifi_type}.{ont_type}.{mapq}.wg_aln-to_T2TXY.ref-cov.tsv',
    output:
        new_names = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.txt',
        name_maps = multiext(
            'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names',
            '.otn-map.tsv', '.otn-map.json', '.otn-map.sed',
            '.nto-map.tsv', '.nto-map.json', '.nto-map.sed'
        )
    wildcard_constraints:
        chrom = '(chrY|chrX)',
        mapq = '(ha|na)'
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
        #wg_assm = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.fasta',
        wg_assm = select_whole_genome_assembly,
        ren_y_json = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.names.otn-map.json',
        ren_x_json = 'output/subset_wg/15_order_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrX.names.otn-map.json',
    output:
        ren_assm = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.fasta',
        sub_y_assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta',
        sub_x_assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrX.fasta',
    wildcard_constraints:
        mapq = '(na|ha)'
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


rule dump_contig_name_mapping_files:
    input:
        fasta = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.fasta'
    output:
        final_names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.txt',
        tsv_otn = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.otn-map.tsv',
        tsv_nto = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.nto-map.tsv',
        json_otn = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.otn-map.json',
        json_nto = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.nto-map.json',
        sed_otn = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.otn-map.sed',
        sed_nto = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.nto-map.sed',
    wildcard_constraints:
        chrom = '(chrY|chrX)',
        sample = SAMPLE_NAME_CONSTRAINT,
        mapq = '(na|ha)'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import json

        new_contig_names = []
        with open(input.fasta, 'r') as fasta:
            for line in fasta:
                if not line.startswith('>'):
                    continue
                new_name = line.strip()[1:]  # ignore leading ">"
                old_name = new_name.split('.')[-2]  # -1 is sample name
                assert new_name.split('.')[-1] == wildcards.sample
                new_contig_names.append((new_name, old_name))

        # dump TSV output
        otn_file = output.tsv_otn
        nto_file = output.tsv_nto
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f'{old_name}\t{new_name}\n')
                    _ = nto.write(f'{new_name}\t{old_name}\n')

        # dump JSON output
        otn_file = output.json_otn
        nto_file = output.json_nto

        otn_map = dict((old_name, new_name) for new_name, old_name in new_contig_names)
        with open(otn_file, 'w') as otn:
            json.dump(otn_map, otn, indent=1)

        nto_map = dict(new_contig_names)
        with open(nto_file, 'w') as nto:
            json.dump(nto_map, nto, indent=1)

        # dump SED output
        otn_file = output.sed_otn
        nto_file = output.sed_nto
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f's/\\b{old_name}\\b/{new_name}/g\n')
                    _ = nto.write(f's/\\b{new_name}\\b/{old_name}/g\n')

        # dump name listing
        with open(output.final_names, 'w') as listing:
            for new_name, _ in new_contig_names:
                _ = listing.write(f'{new_name}\n')
    # END OF RUN BLOCK


rule dump_ref_name_mapping_files:
    """
    Annoying code duplication to process
    chrY reference files for certain steps
    of the pipeline - late incoming request...
    """
    input:
        fasta = 'references_derived/{sample}_{chrom}.fasta'
    output:
        final_names = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.txt',
        tsv_otn = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.otn-map.tsv',
        tsv_nto = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.nto-map.tsv',
        json_otn = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.otn-map.json',
        json_nto = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.nto-map.json',
        sed_otn = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.otn-map.sed',
        sed_nto = 'output/subset_wg/25_name_mappings/{sample}.HIFIRW.ONTUL.na.{chrom}.names.nto-map.sed',
    wildcard_constraints:
        chrom = 'chrY',
        sample = '(T2T|GRCh38)'
    run:
        import json

        new_contig_names = []
        with open(input.fasta, 'r') as fasta:
            for line in fasta:
                if not line.startswith('>'):
                    continue
                new_name = line.strip()[1:]  # ignore leading ">"
                old_name = new_name
                new_contig_names.append((new_name, old_name))

        # dump TSV output
        otn_file = output.tsv_otn
        nto_file = output.tsv_nto
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f'{old_name}\t{new_name}\n')
                    _ = nto.write(f'{new_name}\t{old_name}\n')

        # dump JSON output
        otn_file = output.json_otn
        nto_file = output.json_nto

        otn_map = dict((old_name, new_name) for new_name, old_name in new_contig_names)
        with open(otn_file, 'w') as otn:
            json.dump(otn_map, otn, indent=1)

        nto_map = dict(new_contig_names)
        with open(nto_file, 'w') as nto:
            json.dump(nto_map, nto, indent=1)

        # dump SED output
        otn_file = output.sed_otn
        nto_file = output.sed_nto
        with open(otn_file, 'w') as otn:
            with open(nto_file, 'w') as nto:
                for new_name, old_name in new_contig_names:
                    assert old_name in new_name
                    _ = otn.write(f's/\\b{old_name}\\b/{new_name}/g\n')
                    _ = nto.write(f's/\\b{new_name}\\b/{old_name}/g\n')

        # dump name listing
        with open(output.final_names, 'w') as listing:
            for new_name, _ in new_contig_names:
                _ = listing.write(f'{new_name}\n')
    # END OF RUN BLOCK



rule rename_verkko_coverage_tables:
    input:
        wg_hifi_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.hifi-coverage.csv',
        wg_ont_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.ont-coverage.csv',
        ren_y_sed = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.sed',
        sub_y_bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.bed',
        ren_x_sed = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrX.names.otn-map.sed',
        sub_x_bed = 'output/subset_wg/10_find_contigs/{sample}.{hifi_type}.{ont_type}.na.chrX.bed',
    output:
        ren_hifi_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.hifi-coverage.csv',
        ren_ont_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.ont-coverage.csv',
        merge_sed = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY-chrX.names.otn-map.sed',
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
        paf = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.wg_aln-to_{reference}.paf.gz',
        names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.names.txt'
    output:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'zgrep -w -F -f {input.names} {input.paf} | sort -k1 -k3n,4n | gzip > {output.paf}'


rule cache_contig_alignments:
    """
    Store all primary contig alignments
    for later plotting (more efficient binning)
    """
    input:
        paf = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY_aln-to_{reference}.paf.gz'
    output:
        hdf = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY_aln-to_{reference}.cache.h5'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
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
        df = df.loc[df['tag_tp'].str.lower().isin(['tp:a:p', 'tp:a:i']), :].copy()

        target_size = df.at[df.index[0], 'target_length']
        ctg_cov = np.zeros(target_size, dtype=np.int8)
    
        for tstart, tend in df[['target_start', 'target_end']].itertuples(index=False):
            ctg_cov[tstart:tend] += 1
            if ctg_cov.max() > 254:
                raise
        
        with pd.HDFStore(output.hdf, mode='w', complib='blosc', complevel=9) as hdf:
            hdf.put(wildcards.sample, pd.Series(ctg_cov), format='fixed')
    # END OF RUN BLOCK


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
        names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt'
    output:
        bam = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.na.{chrom}_aln-to_{reference}.bam'
    wildcard_constraints:
        chrom = '(chrX|chrY)'
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
        names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.{chrom}.names.txt'
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


rule cache_read_to_assembly_aln:
    input:
        hifi = 'output/subset_wg/40_extract_rdaln/{sample}.HIFIRW_aln-to_{hifi_type}.{ont_type}.na.{chrom}.paf.gz',
        ont = 'output/subset_wg/40_extract_rdaln/{sample}.ONTUL_aln-to_{hifi_type}.{ont_type}.na.{chrom}.paf.gz',
        fai = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.fai',
    output:
        hdf = 'output/subset_wg/40_extract_rdaln/{sample}.READS_aln-to_{hifi_type}.{ont_type}.na.{chrom}.cache.h5',
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
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

        contigs = []
        with open(input.fai, 'r') as faidx:
            for line in faidx:
                ctg, ctg_size = line.split()[:2]
                contigs.append((ctg, int(ctg_size)))

        hifi_aln = pd.read_csv(input.hifi, sep='\t', header=None, names=PAF_COLUMN_NAMES, usecols=PAF_USE_COLS)
        hifi_aln = hifi_aln.loc[hifi_aln['tag_tp'].str.lower().isin(['tp:a:p', 'tp:a:i']), :].copy()

        ont_aln = pd.read_csv(input.ont, sep='\t', header=None, names=PAF_COLUMN_NAMES, usecols=PAF_USE_COLS)
        ont_aln = ont_aln.loc[ont_aln['tag_tp'].str.lower().isin(['tp:a:p', 'tp:a:i']), :].copy()

        ordered_contigs = []
        for ctg_order, (ctg_name, ctg_size) in enumerate(sorted(contigs), start=1):
            order_key = f'CTG{ctg_order:03}'
            ordered_contigs.append((order_key, ctg_name, ctg_size))

            hifi_cov = np.zeros(ctg_size, dtype=np.int16)
            ont_cov = np.zeros(ctg_size, dtype=np.int16)

            hifi_sub = hifi_aln.loc[hifi_aln['target_chrom'] == ctg_name, ['target_start', 'target_end']]
            for s, e in hifi_sub.itertuples(index=False):
                hifi_cov[s:e] += 1
            
            ont_sub = ont_aln.loc[ont_aln['target_chrom'] == ctg_name, ['target_start', 'target_end']]
            for s, e in ont_sub.itertuples(index=False):
                ont_cov[s:e] += 1
            
            with pd.HDFStore(output.hdf, mode='a', complib='blosc', complevel=9) as hdf:
                # contig names are not "Python identifiers", avoid ugly warnings here
                key = f'{wildcards.sample}/{order_key}/HIFI'
                hdf.put(key, pd.Series(hifi_cov))
                key = f'{wildcards.sample}/{order_key}/ONT'
                hdf.put(key, pd.Series(ont_cov))

        ordered_contigs = pd.DataFrame.from_records(ordered_contigs, columns=['order_key', 'contig_name', 'contig_size'])
        with pd.HDFStore(output.hdf, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put('contigs', ordered_contigs, format='fixed')
    # END OF RUN BLOCK


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
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:59'
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
        names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.names.txt',
        tsv = 'output/motif_search/10_norm/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.{motif}.norm.tsv',
        bed = 'output/motif_search/10_norm/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.{motif}.norm-hiq.bed',
    output:
        tsv = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.{motif}.norm.tsv',
        bed = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.{motif}.norm-hiq.bed',
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
        fasta = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta',
        bed = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.{motif}.norm-hiq.bed',
    output:
        fasta = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.{motif}.hiq-seq.fasta',
    wildcard_constraints:
        sample = SAMPLE_NAME_CONSTRAINT
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'seqtk subseq {input.fasta} {input.bed} > {output.fasta}'


rule extract_motif_hit_sequences_in_reference:
    """
    Just exists for the motif annotation in T2T-Y
    """
    input:
        fasta = 'references_derived/{sample}_chrY.fasta',
        bed = 'output/motif_search/10_norm/20_refseq/{sample}.HIFIRW.ONTUL.na.chrY.{motif}.norm-hiq.bed',
    output:
        fasta = 'output/motif_search/10_norm/20_refseq/{sample}.HIFIRW.ONTUL.na.chrY.{motif}.hiq-seq.fasta',
    wildcard_constraints:
        sample = 'T2T'
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


rule cache_read_to_reference_aln:
    input:
        hifi = 'output/subset_wg/60_subset_rdref/{sample}.HIFIRW_aln-to_{reference}.{chrom}.paf.gz',
        ont = 'output/subset_wg/60_subset_rdref/{sample}.ONTUL_aln-to_{reference}.{chrom}.paf.gz',
        fai = 'references_derived/T2T_{chrom}.fasta.fai',
    output:
        hdf = 'output/subset_wg/60_subset_rdref/{sample}.READS_aln-to_{reference}.{chrom}.cache.h5',
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
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

        contigs = []
        with open(input.fai, 'r') as faidx:
            for line in faidx:
                ctg, ctg_size = line.split()[:2]
                contigs.append((ctg, int(ctg_size)))

        hifi_aln = pd.read_csv(input.hifi, sep='\t', header=None, names=PAF_COLUMN_NAMES, usecols=PAF_USE_COLS)
        hifi_aln = hifi_aln.loc[hifi_aln['tag_tp'].str.lower().isin(['tp:a:p', 'tp:a:i']), :].copy()

        ont_aln = pd.read_csv(input.ont, sep='\t', header=None, names=PAF_COLUMN_NAMES, usecols=PAF_USE_COLS)
        ont_aln = ont_aln.loc[ont_aln['tag_tp'].str.lower().isin(['tp:a:p', 'tp:a:i']), :].copy()

        ordered_contigs = []
        for ctg_order, (ctg_name, ctg_size) in enumerate(sorted(contigs), start=1):
            order_key = f'CTG{ctg_order:03}'
            ordered_contigs.append((order_key, ctg_name, ctg_size))

            hifi_cov = np.zeros(ctg_size, dtype=np.int16)
            ont_cov = np.zeros(ctg_size, dtype=np.int16)

            hifi_sub = hifi_aln.loc[hifi_aln['target_chrom'] == ctg_name, ['target_start', 'target_end']]
            for s, e in hifi_sub.itertuples(index=False):
                hifi_cov[s:e] += 1
            
            ont_sub = ont_aln.loc[ont_aln['target_chrom'] == ctg_name, ['target_start', 'target_end']]
            for s, e in ont_sub.itertuples(index=False):
                ont_cov[s:e] += 1
            
            with pd.HDFStore(output.hdf, mode='a', complib='blosc', complevel=9) as hdf:
                # contig names are not "Python identifiers", avoid ugly warnings here
                key = f'{wildcards.sample}/{order_key}/HIFI'
                hdf.put(key, pd.Series(hifi_cov))
                key = f'{wildcards.sample}/{order_key}/ONT'
                hdf.put(key, pd.Series(ont_cov))

        ordered_contigs = pd.DataFrame.from_records(ordered_contigs, columns=['order_key', 'contig_name', 'contig_size'])
        with pd.HDFStore(output.hdf, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put('contigs', ordered_contigs, format='fixed')
    # END OF RUN BLOCK