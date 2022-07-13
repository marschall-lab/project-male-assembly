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
            ('autosomal', 'haploid', autosomes, num_auto),
            ('autosomal', 'diploid', autosomes * 2, num_auto * 2),
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


rule norm_expert_seqclasses_file:
    """
    This rule is needed to first get rid
    of the CR line terminators from Pille's
    seq. classes annotation files.

    cut: sometime strand info is present
    """
    input:
        bed_dos = ancient(
            pl.Path(f"{config['path_root_share_working']}",
            "assemblies/verkko_1.0_release/chrY/{sample}.HIFIRW.ONTUL.na.chrY_SeqClasses.bed")
        )
    output:
        bed_unix = 'references_derived/{sample}.HIFIRW.ONTUL.na.chrY.seqclasses.bed'
    shell:
        'cat {input.bed_dos} | dos2unix | cut -f 1-4 > {output.bed_unix}'


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
                    checksums.append((partnum, filename, md5sum))
            df = pd.DataFrame.from_records(checksums, columns=['part', 'filename', 'md5'])
            hdf.put('checksums', df, format='fixed')
    # END OF RUN BLOCK
