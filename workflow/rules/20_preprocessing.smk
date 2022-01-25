import pathlib as pl


rule merge_hg002_chry_draft:
    """
    Add current T2T/chrY draft to
    T2T/chm13 assembly
    """
    input:
        genome = '/gpfs/project/projects/medbioinf/data/references/T2Tv11_T2TC_chm13.fasta',
        chry = '/gpfs/project/projects/medbioinf/data/references/hg002.chrY.v2.fasta'
    output:
        genome = '/gpfs/project/projects/medbioinf/data/references/T2Tv11_hg002Yv2_chm13.fasta'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import shutil
        import io

        shutil.copy(input.genome, output.genome)

        out_buffer = io.StringIO()
        _ = out_buffer.write('>chrY\n')
        with open(input.chry, 'r') as fasta_in:
            _ = fasta_in.readline()
            _ = out_buffer.write(fasta_in.read())

        with open(output.genome, 'a') as fasta_out:
            _ = fasta_out.write(out_buffer.getvalue())


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
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.{seq_type}.gz',
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
        'seqtk subseq {input.reads} {input.names} | pigz -p {threads} --best > {output}'


rule merge_sex_chrom_reads:
    input:
        chrx = 'output/read_subsets/chrX/{sample_info}_{sample}_{read_type}.chrX-reads.{mapq}.{seq_type}.gz',
        chry = 'output/read_subsets/chrY/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.{seq_type}.gz',
    output:
        'output/read_subsets/chrXY/{sample_info}_{sample}_{read_type}.chrXY-reads.{mapq}.{seq_type}.gz',
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
        hifiec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_HIFIEC.{chrom}-reads.{mapq}.{seq_type}.gz',
        ontec = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTEC.{chrom}-reads.{mapq}.{seq_type}.gz',
    output:
        'output/read_subsets/{chrom}/{sample_info}_{sample}_OHEC.{chrom}-reads.{mapq}.{seq_type}.gz',
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
        'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.{seq_type}.gz',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample = '(NA193N7|NA193NN|AFR4MIX)'
    threads: config['num_cpu_low']
    resources:
        walltime = lambda wildcards, attempt: f'{attempt * attempt:02}:00:00',
    shell:
        'pigz -d -c {input.reads} | pigz --best -p {threads} > {output}'
