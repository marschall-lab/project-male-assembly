
# For files that are not dependant on a Verkko version,
# do not copy files again if the version on the Globus
# share has the same size as the source file. At the time
# of writing, this is only relevant for reference files
# and read sets.
SKIP_SIZE = config.get('populate_share_skip_size', False)

# TODO
# entire module should be switched to rsync for simplicity

rule copy_references:
    input:
        grch38_noalt = 'references_derived/GRCh38_noalt.fasta',
        grch38_noalt_idx = 'references_derived/GRCh38_noalt.fasta.fai',
        grch38_chry = 'references_derived/GRCh38_chrY.fasta',
        t2t_chm13 = 'references_derived/T2T_chm13_122XM.fasta',
        t2t_chm13_idx = 'references_derived/T2T_chm13_122XM.fasta.fai',
        t2t_hg002 = 'references_derived/T2T_122XYM.fasta',
        t2t_hg002_idx = 'references_derived/T2T_122XYM.fasta.fai',
        t2t_chry = 'references_derived/T2T_chrY.fasta',
        motifs = expand(
            'references_derived/{motif}.fasta',
            motif=config['motif_search']
        )
    output:
        ok = 'output/share/references.ok'
    run:
        import os
        import pathlib as pl
        import shutil as sh
    
        share_path = pl.Path(config['path_root_share_references']).resolve()

        check_file = ''
        for ref_file in input:
            make_copy = True
            source = pl.Path(ref_file)
            dest = share_path / source.name
            if dest.is_file() and SKIP_SIZE:
                source_size = os.stat(source).st_size
                dest_size = os.stat(dest).st_size
                make_copy = source_size != dest_size
            if make_copy:
                sh.copy2(source, dest)
            check_file += f'{source}\t{dest}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_chromosome_readsets:
    input:
        fasta = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.fasta.gz'
    output:
        ok = 'output/share/read_subsets/{sample_info}_{sample}_{read_type}.{chrom}-reads.{mapq}.copied.ok'
    run:
        import os
        import pathlib as pl
        import shutil as sh

        share_path = pl.Path(config['path_root_share_working']).resolve()
        readsets_subfolder = share_path / pl.Path(f'readsets/{wildcards.chrom}/{wildcards.read_type}')
        readsets_subfolder.mkdir(parents=True, exist_ok=True)

        source = pl.Path(input.fasta)
        destination = readsets_subfolder / source.name
        if destination.is_file() and SKIP_SIZE:
            source_size = os.stat(source).st_size
            dest_size = os.stat(destination).st_size
            if source_size != dest_size:
                sh.copy(source, destination)
        else:
            sh.copy(source, destination)

        with open(output.ok, 'w') as dump:
            _ = dump.write(f'{source}\t{destination}\n')
    # END OF RUN BLOCK


rule copy_verkko_assemblies:
    """
    Collect additional files on-the-fly
    HPC outputs are not useful for downstream analyses (at the moment)
    """
    input:
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.info',
        linear = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
    output:
        ok = 'output/share/assemblies/verkko_{major}_{minor}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    run:
        import pathlib as pl
        import shutil as sh

        skip_hpc_outputs = True

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'assemblies',
            wildcards.chrom
        )

        check_file = ''

        assembly_pattern = pl.Path(input.linear).with_suffix('').name
        assembly_files = pl.Path(input.linear).parent.glob(f'./{assembly_pattern}*')

        for assm_file in assembly_files:
            new_name = pl.Path(assembly_id + '.' + assm_file.name)
            if 'compressed' in str(new_name) and skip_hpc_outputs:
                continue
            dest_file = verkko_subfolder / new_name
            sh.copy2(assm_file, dest_file)
            check_file += f'{assm_file}\t{dest_file}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_contig_to_ref_alignments:
    input:
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.info',
        bam = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.sort.bam',
        bai = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.sort.bam.bai',
        paf = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz',
        bed = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed',
        cov = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.ref-cov.tsv',
    output:
        ok = 'output/share/alignments/contigs-to-ref/verkko_{major}_{minor}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.copied.ok'
    run:
        import pathlib as pl
        import shutil as sh

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/contigs-to-ref',
            wildcards.chrom
        )
        check_file = ''

        for source in [input.bam, input.bai, input.paf, input.bed, input.cov]:
            source_path = pl.Path(source)
            dest_path = verkko_subfolder / source_path.name
            sh.copy2(source_path, dest_path)
            check_file += f'{source_path}\t{dest_path}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_reads_to_assm_alignments:
    input:
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.info',
        bam = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.bam',
        bai = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.bam.bai',
        paf = 'output/alignments/reads-to-assm/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.paf.gz'
    output:
        ok = 'output/share/alignments/reads-to-assm/verkko_{major}_{minor}/{sample_info}_{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    run:
        import pathlib as pl
        import shutil as sh

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/reads-to-assm',
            wildcards.chrom
        )
        check_file = ''

        for source in [input.bam, input.bai, input.paf]:
            source_path = pl.Path(source)
            dest_path = verkko_subfolder / source_path.name
            sh.copy2(source_path, dest_path)
            check_file += f'{source_path}\t{dest_path}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK
