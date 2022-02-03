import pathlib as pl


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
        'grep -v decoy {input.index} > {output.listing}'


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
        'seqtk subseq {input.genome} {input.chrom_list} > {output.genome}'


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
        'echo chrY > tmp.name && '
        'seqtk subseq {input.genome} tmp.name > {output.chrom} && '
        'rm tmp.name'



###################################################################
# Below: extract subsets of reads aligning to specific chromosomes
# here: T2T/chrX and T2T/chrY are relevant
# NB: the whole-genome alignment step including the generation of
# the output cache is implemented as part of the HGSVC ONT-QC
# pipeline; should be moved here at some point
###################################################################


rule extract_aligned_chry_read_names:
    input:
        cached_aln = select_alignment_cache_file
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.mq60.txt',
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.mq00.txt'
    run:
        import pandas as pd
        with pd.HDFStore(input.cached_aln, 'r') as hdf:
            chrom_aln = hdf[wildcards.chrom]
            highq_reads = chrom_aln.loc[chrom_aln['mapq'] == 60, 'read_name'].unique().tolist()
            with open(output[0], 'w') as dump:
                _ = dump.write('\n'.join(sorted(highq_reads)) + '\n')
            all_reads = chrom_aln['read_name'].unique().tolist()
            with open(output[1], 'w') as dump:
                _ = dump.write('\n'.join(sorted(all_reads)) + '\n')


rule extract_aligned_chrom_read_sequences:
    input:
        names = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.txt',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type],
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.fasta.gz',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample = CONSTRAINT_REGULAR_SAMPLES,
        read_type = '(HIFIEC|HIFIAF|ONTUL|ONTEC)',
        chrom = '(chrX|chrY)'
    threads: config['num_cpu_low']
    resources:
        walltime = lambda wildcards, attempt: f'{4 * attempt:02}:00:00',
    shell:
        'seqtk subseq {input.reads} {input.names} | seqtk seq -A -C | pigz -p {threads} --best > {output}'


rule merge_sex_chrom_reads:
    input:
        chrx = 'output/read_subsets/chrX/{sample_info}_{sample}_{read_type}.chrX-reads.{mapq}.fasta.gz',
        chry = 'output/read_subsets/chrY/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.fasta.gz',
    output:
        'output/read_subsets/chrXY/{sample_info}_{sample}_{read_type}.chrXY-reads.{mapq}.fasta.gz',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF|ONTUL|ONTEC)',
        sample = CONSTRAINT_REGULAR_SAMPLES,
    threads: config['num_cpu_low']
    resources:
        walltime = lambda wildcards, attempt: f'{attempt:02}:00:00',
    shell:
        'pigz -d -c {input.chrx} {input.chry} | pigz --best -p {threads} > {output}'


rule merge_read_types:
    input:
        hifiec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_HIFIEC.{chrom}-reads.{mapq}.fasta.gz',
        ontec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTEC.{chrom}-reads.{mapq}.fasta.gz',
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_OHEC.{chrom}-reads.{mapq}.fasta.gz',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        chrom = '(chrY|chrX|chrXY)',
        sample = CONSTRAINT_REGULAR_SAMPLES,
    threads: config['num_cpu_low']
    resources:
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:00:00',
    shell:
        'pigz -d -c {input.hifiec} {input.ontec} | pigz --best -p {threads} > {output}'


rule merge_afr_mix_subsets:
    """
    Dedicated rule to merge read sets for the AFR mix sample: NA19317 + NA19347 = NA193N7
    """
    input:
        reads = select_afr_mix_subsets
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.fasta.gz',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample = '(NA193N7|NA193NN|AFR4MIX)'
    threads: config['num_cpu_low']
    resources:
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:00:00',
    shell:
        'pigz -d -c {input.reads} | pigz --best -p {threads} > {output}'
