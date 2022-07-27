import pathlib as pl

localrules: deselect_decoy_chroms_grch38, \
            norm_expert_seqclasses_file, \
            prep_t2t_seq_class_cache_file, \
            dump_reference_genome_sizes

rule t2t_convert_to_ucsc_ids:
    input:
        genbank = 'references/GCA_009914755.3_CHM13_T2T_v1.1.fasta',
        mapping = expand('references/{chr2acc}',
            chr2acc=[
                'GCA_009914755.3_CHM13_T2T_v1.1.primary-chr2acc.tsv',
                'GCA_009914755.3_CHM13_T2T_v1.1.mito-chr2acc.tsv'
            ]
        )
    output:
        local = 'references_derived/T2T_chm13_122XM.fasta'
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt
    run:
        import pandas as pd
        import io

        lut = []
        for mapping in input.mapping:
            df = pd.read_csv(mapping, sep='\t', header=None, names=['short', 'acc'], comment='#')
            lut.append(df)
        lut = pd.concat(lut)
        lut['short'] = lut['short'].apply(lambda x: 'chrM' if x == 'MT' else f'chr{x}')
        lut = {acc:chrom for chrom, acc in lut.itertuples(index=None)}

        out_buffer = io.StringIO()
        with open(input.genbank, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    new_names = [v for k, v in lut.items() if k in line]
                    if len(new_names) != 1:
                        raise ValueError(f'No unique name mapping for input chromosome {line.strip()}')
                    this_chrom = f'>{new_names[0]}\n'
                    out_buffer.write(this_chrom)
                    continue
                out_buffer.write(line)

        with open(output.local, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())
    # END OF RUN BLOCK


rule t2t_add_hg002_to_chm13:
    input:
        genome = 'references_derived/T2T_chm13_122XM.fasta',
        chry = 'references/CP086569.2.fasta'
    output:
        genome = 'references_derived/T2T_122XYM.fasta',
        chry = 'references_derived/T2T_chrY.fasta'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import shutil
        import io

        shutil.copy(input.genome, output.genome)

        out_buffer = io.StringIO()
        _ = out_buffer.write('>chrY\n')
        with open(input.chry, 'r') as fasta_in:
            header = fasta_in.readline()
            if not 'CP086569.2' in header:
                raise ValueError(f'This is not the correct version of chrY: {header.strip()}')
            _ = out_buffer.write(fasta_in.read())

        with open(output.genome, 'a') as fasta_out:
            _ = fasta_out.write(out_buffer.getvalue())

        with open(output.chry, 'w') as fasta_out:
            _ = fasta_out.write(out_buffer.getvalue())
    # END OF RUN BLOCK


rule deselect_decoy_chroms_grch38:
    input:
        index = 'references/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai'
    output:
        listing = 'references_derived/GRCh38_keep_chromosomes.txt'
    shell:
        'grep -v decoy {input.index} | cut -f 1 > {output.listing}'


rule extract_non_decoy_chroms_grch38:
    input:
        genome = 'references/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta',
        chrom_list = 'references_derived/GRCh38_keep_chromosomes.txt'
    output:
        genome = 'references_derived/GRCh38_noalt.fasta'
    conda:
        '../envs/biotools.yaml'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt:02}:59:00',
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
    shell:
        'seqtk subseq {input.genome} {input.chrom_list} | seqtk seq -C > {output.genome}'


rule extract_grch38_chrom_y:
    input:
        genome = 'references_derived/GRCh38_noalt.fasta'
    output:
        chrom = 'references_derived/GRCh38_chrY.fasta'
    conda:
        '../envs/biotools.yaml'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt:02}:59:00',
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
    shell:
        'echo chrY > tmp.name && echo chrY_KI270740v1_random >> tmp.name && '
        'seqtk subseq {input.genome} tmp.name > {output.chrom} && '
        'rm tmp.name'


rule create_reference_bed_file:
    """
    This rule adds the chrY sequence classes (PAR, XDR etc.) to the main
    chromosomes in the resulting BED file. This can be used to summarize
    coverages per region etc.
    """
    input:
        fai = lambda wildcards: select_reference_genome(wildcards.reference, True),
        seq_classes = lambda wildcards: f'references_derived/{config["reference_y_seq_classes"][wildcards.reference]}.bed'
    output:
        bed = 'references_derived/{reference}.bed'
    wildcard_constraints:
        reference = '(' + '|'.join(sorted(config['reference_genomes'].keys())) + ')'
    run:
        bed_out = []
        with open(input.fai, 'r') as fasta_index:
            for line in fasta_index:
                chrom, chrom_size = line.strip().split()[:2]
                if chrom == 'chrY':
                    # note here: GRCh38 also contains an
                    # unplaced chrY contig, that is kept
                    # as-is b/c no sequence classes exist
                    # for that one
                    continue
                bed_out.append((chrom, '0', chrom_size, chrom))

        with open(input.seq_classes) as table:
            for line in table:
                chrom, start, end, name = line.strip().split()[:4]
                bed_out.append((chrom, start, end, name))

        with open(output.bed, 'w') as dump:
            _ = dump.write('\n'.join(['\t'.join(record) for record in bed_out]) + '\n')
    # END OF RUN BLOCK


rule prep_t2t_seq_class_cache_file:
    """
    Special rule to prepare the T2T-Y sequence class
    annotations to be suitable for plotting
    """
    input:
        bed = 'references_derived/T2T.chrY-seq-classes.bed'
    output:
        tsv = 'references_derived/T2T.chrY-seq-classes.tsv'
    run:
        import matplotlib.colors as mcol
        import pandas as pd

        names = ['chrom', 'start', 'end', 'name', 'hex_color', 'full_name']
        df = pd.read_csv(input.bed, sep='\t', header=None, names=names)
        # next line: manual fix for OBO
        df.loc[df['name'] == 'PAR2', 'end'] -= 1
        assert df.at[df.index[-1], 'end'] == 62460029
        # convert HEX to RGBA
        df[['red', 'green', 'blue', 'alpha']] = 0.
        rgb_colors = []
        for row in df.itertuples(index=False):
            r,g,b,a = mcol.to_rgba(row.hex_color)
            rgb_colors.append([r, g, b, a])
        df[['red', 'green', 'blue', 'alpha']] = rgb_colors
    
        df.to_csv(output.tsv, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule dump_reference_genome_sizes:
    input:
        fai = lambda wildcards: select_reference_genome(wildcards.reference, True),
    output:
        table = 'references_derived/{reference}.gsize.tsv'
    run:
        import pandas as pd

        autosomes = 0
        num_auto = 0
        sex_male = 0
        num_male = 0
        sex_female = 0
        num_female = 0

        with open(input.fai, 'r') as faidx:
            for line in faidx:
                chrom, chrom_size = line.split()[:2]
                if 'chrY' in chrom:
                    sex_male += int(chrom_size)
                    num_male += 1
                elif 'chrX' in chrom:
                    sex_female += int(chrom_size)
                    num_female += 1
                else:
                    autosomes += int(chrom_size)
                    num_auto += 1

        gsizes = [
            ('auto_mito', 'haploid', autosomes, num_auto),
            ('auto_mito', 'diploid', autosomes * 2, num_auto * 2),
            ('female', 'haploid', autosomes + sex_female, num_auto + num_female),
            ('female', 'diploid', (autosomes + sex_female) * 2, (num_auto + num_female) * 2),
            ('male', 'diploid', autosomes * 2 + sex_male + sex_female, num_auto * 2 + num_male + num_female),
            ('unisex', 'linear', autosomes + sex_male + sex_female, num_auto + num_male + num_female)
        ]
        gsizes = pd.DataFrame.from_records(gsizes, columns=['karyotype', 'ploidy', 'size', 'chromosomes'])
        gsizes.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


#######################################################################
# Special rule:
# - sync manually created motif and annotation files from the Globus
# share "references/" to the local reference folder to make rules
# less bloated with path lookups etc. and to avoid that a Globus mess
# triggers a pipeline rerun
# - this mainly applies to reference files created via expert curation
#######################################################################


rule sync_expert_reference_file:
    input:
        infile = ancient(f'{config["path_root_share_references"]}/' + '{file_name}')
    output:
        outfile = 'references_derived/{file_name}'
    wildcard_constraints:
        file_name = '(' + '|'.join(config['expert_annotations']) + ')'
    shell:
        'rsync --checksum {input.infile} {output.outfile}'
            ' && '
        'md5sum {output.outfile} > {output.outfile}.md5'


rule normalize_expert_seqclasses_file:
    """
    This rule is needed to first get rid
    of the CR line terminators from Pille's
    seq. classes annotation files.

    cut: sometime strand info is present

    2022-07-27
    adapt this rule to manual preprocessing of the "seqclasses"
    files to fix the following issues:
    - check for non-existing sequence class names
    - add indicator columns to make subsequent data aggregation easier
    - Pandas-internal: input files may contain CR line terminators
    """
    input:
        t2t = 'references_derived/T2T.chrY-seq-classes.tsv',
        bed_dos = ancient(
            pl.Path(f"{config['path_root_share_working']}",
            "assemblies/verkko_1.0_release/chrY/{sample}.HIFIRW.ONTUL.na.chrY_SeqClasses.bed")
        )
    output:
        bed_generic = 'references_derived/seqclasses/{sample}.HIFIRW.ONTUL.na.chrY.generic-seqcls.bed',
        bed_specific = 'references_derived/seqclasses/{sample}.HIFIRW.ONTUL.na.chrY.specific-seqcls.bed',
        table_complete = 'references_derived/seqclasses/{sample}.HIFIRW.ONTUL.na.chrY.seqclasses.tsv',
    run:
        import pandas as pd
        import re as re

        t2t = pd.read_csv(input.t2t, sep='\t', header=0)
        t2t['t2t_length'] = t2t['end'] - t2t['start']
        region_indices = dict((row.name, row.Index) for row in t2t.itertuples())

        assm = pd.read_csv(
            input.bed_dos, sep='\t', header=None,
            names=['contig', 'start', 'end', 'name'],
            usecols=[0, 1, 2, 3]  # sometimes, strand is included
        )
        assm['length'] = assm['end'] - assm['start']

        # indicator for unplaced contigs
        match_unplaced = re.compile('unplaced')
        assm['is_unplaced'] = 0
        assm['is_unplaced'] = (assm['name'].str.contains(match_unplaced)).astype(int)

        # indicator for split annotations, also triggers for enum unplaced...
        match_split_num = re.compile('_[0-9]+$')
        assm['is_split'] = 0
        assm['is_split'] = (assm['name'].str.contains(match_split_num)).astype(int)
        # ...correct for that
        assm.loc[assm['is_unplaced'] == 1, 'is_split'] = 0

        # determine contig span
        assm[['start_seqclass', 'end_seqclass']] = ['n/a', 'n/a']
        assm['start_seqclass'] = assm['contig'].apply(lambda x: x.split('.')[3].split('-')[0])
        assm['end_seqclass'] = assm['contig'].apply(lambda x: x.split('.')[3].split('-')[1])
        assm['start_idx'] = assm['start_seqclass'].replace(region_indices)
        assm['end_idx'] = assm['end_seqclass'].replace(region_indices)

        # heuristic to check for unidentifiable names:
        # after ignoring enumerated seq. classes or unplaced contigs,
        # all that's left should have a proper name
        select_not_split = assm['is_split'] == 0
        select_is_placed = assm['is_unplaced'] == 0
        normal_names = assm.loc[select_not_split & select_is_placed, 'name'].values
        error_names = set(normal_names) - set(t2t['name'])
        if error_names:
            raise ValueError(error_names)

        match_seqclass = '^(' + '|'.join(sorted(t2t['name'])) + ')'
        assm['seqclass'] = assm['name'].str.extract(match_seqclass)
        assm['seqclass_idx'] = assm['seqclass'].replace(region_indices)

        # aggregate contiguity information per sequence class; makes
        # many operations downstream trivial
        contiguous_regions = []
        for seqclass, annotations in assm.groupby('seqclass'):
            num_contigs = annotations['contig'].nunique()
            assm_length = annotations['length'].sum()
            ref_length = t2t.loc[t2t['name'] == seqclass, 't2t_length'].values[0]
            assm_length_pct = round(assm_length / ref_length * 100, 2)
            is_contiguous = 1
            if (annotations['is_unplaced'] == 1).any():
                is_contiguous = 0
            if (annotations['is_split'] == 1).any() and num_contigs > 1:
                is_contiguous = 0
            if num_contigs == 1:
                left_idx = annotations['start_idx'].values[0]
                right_idx = annotations['end_idx'].values[0]
                my_idx = annotations['seqclass_idx'].values[0]
                if left_idx < my_idx < right_idx:
                    is_contiguous = 1
                elif seqclass in ['PAR1', 'PAR2']:
                    if assm_length_pct > 95:
                        is_contiguous = 1
                else:
                    is_contiguous = 0
            contiguous_regions.append(
                (
                    seqclass,
                    num_contigs,
                    assm_length,
                    ref_length,
                    assm_length_pct,
                    is_contiguous
                )
            )
            
        cr_columns = [
            'seqclass', 'assm_contigs_num',
            'assm_length_bp', 'ref_length_bp',
            'assm_length_pct', 'is_contiguous'
        ]
        contiguous_regions = pd.DataFrame.from_records(
            contiguous_regions,
            columns=cr_columns
        )

        assm = assm.merge(contiguous_regions, left_on='seqclass', right_on='seqclass', how='outer')
        assert pd.notnull(assm).all(axis=0).all()
        assm.sort_values(['contig', 'start', 'end'], ascending=True, inplace=True)

        assm.to_csv(output.table_complete, sep='\t', header=True, index=False, line_terminator="\n")
        assm[['contig', 'start', 'end', 'seqclass']].to_csv(
            output.bed_generic, sep='\t', header=False, index=False, line_terminator="\n"
        )
        assm[['contig', 'start', 'end', 'name']].to_csv(
            output.bed_specific, sep='\t', header=False, index=False, line_terminator="\n"
        )
    # END OF RUN BLOCK


rule compute_read_stats:
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample][wildcards.reads][int(wildcards.partnum)]
    output:
        table = 'output/stats/reads/parts/{sample}.{reads}.part{partnum}.tsv.gz'
    conda:
        '../envs/biotools.yaml'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*6:02}:59:59'
    shell:
        'seqtk comp {input.reads} | cut -f 1-6 | pigz -p {threads} --best > {output.table}'


rule compute_read_checksum:
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample][wildcards.reads][int(wildcards.partnum)]
    output:
        md5 = 'output/checksums/{sample}.{reads}.part{partnum}.md5'
    conda:
        '../envs/biotools.yaml'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*6:02}:59:59'
    shell:
        'md5sum {input.reads} > {output.md5}'


rule merge_read_stats:
    input:
        tables = lambda wildcards: expand(
            'output/stats/reads/parts/{{sample}}.{{reads}}.part{partnum}.tsv.gz',
            partnum=list(range(0,len(SAMPLE_DATA[wildcards.sample][wildcards.reads])))
        ),
        md5 = lambda wildcards: expand(
            'output/checksums/{{sample}}.{{reads}}.part{partnum}.md5',
            partnum=list(range(0,len(SAMPLE_DATA[wildcards.sample][wildcards.reads])))
        )
    output:
        hdf = 'output/stats/reads/cached/{sample}.{reads}.read-stats.h5'
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:59'
    run:
        import pandas as pd
        import pathlib as pl

        checksums = []
        stat_cols = ['read_name', 'read_length', 'num_A', 'num_C', 'num_G', 'num_T']
        with pd.HDFStore(output.hdf, mode='w', complib='blosc', complevel=9) as hdf:
            for stats_file, chk_file in zip(input.tables, input.md5):
                partnum = chk_file.rsplit('.', 2)[-2]
                assert partnum.startswith('part')
                assert partnum in stats_file
                df = pd.read_csv(stats_file, sep='\t', header=None, names=stat_cols)
                hdf.put(f'stats/{partnum}', df, format='fixed')
                with open(chk_file, 'r') as text:
                    md5sum, filename = text.read().strip().split()
                checksums.append((partnum, pl.Path(filename).name, md5sum))
            df = pd.DataFrame.from_records(checksums, columns=['part', 'filename', 'md5'])
            hdf.put('checksums', df, format='fixed')
    # END OF RUN BLOCK


rule compute_read_coverage_stats:
    input:
        gsize = 'references_derived/T2TXY.gsize.tsv',
        read_cache = 'output/stats/reads/cached/{sample}.{reads}.read-stats.h5'
    output:
        table = 'output/stats/reads/{sample}.{reads}.read-stats.tsv'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import numpy as np
        import collections as col
        use_sizes = [(('unisex', 'linear'), 'T2TXYM_linear'), (('male', 'diploid'), 'T2TXYM_diploid')]
        gsizes = pd.read_csv(input.gsize, sep="\t", header=0)

        ref_sizes = {}
        for (karyo, ploidy), size_label in use_sizes:
            select_k = gsizes['karyotype'] == karyo
            select_p = gsizes['ploidy'] == ploidy
            selector = select_k & select_p
            ref_size = gsizes.loc[selector, 'size'].values[0]
            ref_sizes[size_label] = ref_size


        reads = []
        with pd.HDFStore(input.read_cache, 'r') as hdf:
            for k in hdf.keys():
                if 'checksums' in k:
                    continue
                read_lengths = hdf[k]['read_length'].values
                reads.append(read_lengths)
        reads = np.concatenate(reads, dtype=np.int32)
        reads.sort()

        read_length_median = reads[reads.size // 2]
        read_length_mean = int(round(reads.mean(), 0))

        num_bp = reads.sum()
        read_length_n50 = reads[reads.cumsum() >= (num_bp // 2)].min()

        stats = col.OrderedDict({
            'sample': wildcards.sample,
            'read_type': wildcards.reads,
            'read_length_N50_bp': read_length_n50,
            'read_length_N50_kbp': int(round(read_length_n50 / 1e3, 0))
        })

        thresholds = [(0, 'geq_0bp'), (15000, 'geq_15kbp'), (1e5, 'geq_100kbp'), (1e6, 'geq_1Mbp')]

        for t, t_label in thresholds:
            num_reads = int((reads >= t).sum())
            stats[f'num_reads_{t_label}'] = num_reads
            if t_label != 'geq_0bp':
                pct_reads = round(num_reads / stats['num_reads_geq_0bp'] * 100, 1)
                stats[f'pct_reads_{t_label}'] = pct_reads
            t_num_bp = reads[reads >= t].sum()
            stats[f'num_bp_{t_label}'] = t_num_bp
            stats[f'num_Gbp_{t_label}'] = round(t_num_bp / 1e9, 2)

            for ref_label, ref_size in ref_sizes.items():
                label = f'cov_{t_label}_{ref_label}'
                cov = int(round(t_num_bp / ref_size, 0))
                stats[label] = cov
        
        stats.update(col.OrderedDict({
            'read_length_min_bp': reads.min(),
            'read_length_median_bp': read_length_median,
            'read_length_median_kbp': int(round(read_length_median / 1e3, 0)),
            'read_length_mean_bp': read_length_mean,
            'read_length_mean_kbp': int(round(read_length_mean / 1e3, 0)),
            'read_length_max_bp': reads.max(),
            'read_length_max_kbp': int(round(reads.max() / 1e3, 0))
        }))

        for ref_label, ref_size in ref_sizes.items():
            stats[f'{ref_label}_bp'] = ref_size
            stats[f'{ref_label}_Gbp'] = round(ref_size / 1e9, 2)

        df = pd.DataFrame.from_records([stats])
        df.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_all_read_stats

rule merge_all_read_stats:
    input:
        tables = expand(
            'output/stats/reads/{sample}.{reads}.read-stats.tsv',
            sample=SAMPLE_NAMES,
            reads=['HIFIRW', 'ONTUL']
        )
    output:
        table = 'output/stats/reads/SAMPLES.READS.read-stats.tsv'
    run:
        import pandas as pd

        merged = []
        for table in input.tables:
            df = pd.read_csv(table, sep='\t', header=0)
            merged.append(df)

        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(['sample', 'read_type'], ascending=True, inplace=True)

        merged.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK
