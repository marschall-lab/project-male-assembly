import pathlib as pl

localrules: deselect_decoy_chroms_grch38

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
        # convert HEX to RGBA
        df[['red', 'green', 'blue', 'alpha']] = 0.
        rgb_colors = []
        for row in df.itertuples(index=False):
            r,g,b,a = mcol.to_rgba(row.hex_color)
            rgb_colors.append([r, g, b, a])
        df[['red', 'green', 'blue', 'alpha']] = rgb_colors
    
        df.to_csv(output.tsv, sep='\t', header=True, index=False)
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
