

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
    output:
        ok = 'output/share/references.ok'
    run:
        import pathlib as pl
        import shutil as sh
    
        share_path = pl.Path(config['path_root_share_references']).resolve()

        check_file = ''
        for ref_file in input:
            source = pl.Path(ref_file)
            dest = share_path / source.name
            sh.copy(source, dest)
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
        import pathlib as pl
        import shutil as sh

        share_path = pl.Path(config['path_root_share_working']).resolve()
        readsets_subfolder = share_path / pl.Path(f'readsets/{wildcards.chrom}/{wildcards.read_type}')
        readsets_subfolder.mkdir(parents=True, exist_ok=True)

        source = pl.Path(input.fasta)
        destination = readsets_subfolder / source.name
        sh.copy(source, destination)

        with open(output.ok, 'w') as dump:
            _ = dump.write(f'{source}\t{destination}\n')
    # END OF RUN BLOCK


rule copy_chromosome_assemblies:
    input:
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verrko.info',
        linear = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        graph = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.gfa',
        stats = multiext(
            'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly',
            '.hifi-coverage.csv', '.layout', '.ont-coverage.csv'
        )
    output:
        ok = 'output/share/assemblies/verkko_{release}_{commit}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    run:
        import json
        import pathlib as pl
        import shutil as sh

        with open(input.version, 'r') as js_dump:
            verkko = json.load(js_dump)
            release = verkko['verkko_release']
            commit = verkko['verkko_commit']

        share_path = pl.Path(config['path_root_share_working']).resolve()
        verkko_subfolder = share_path / pl.Path(f'assemblies/verkko_{release}_{commit}/{wildcards.chrom}/{wildcards.hifi_type}')
        verkko_subfolder.mkdir(parents=True, exist_ok=True)

        check_file = ''

        for source in [input.linear, input.graph] + list(input.stats):
            # not sure if easier to separate by output type...
            new_name = pl.Path(source.replace('/assembly', '.assembly')).name
            source_path = pl.Path(source)
            dest_path = verkko_subfolder / new_name
            sh.copy(source_path, dest_path)
            check_file += f'{source_path}\t{dest_path}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_chromosome_contig_alignments:
    input:
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verrko.info',
        bam = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.sort.bam',
        bai = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.sort.bam.bai',
        paf = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    output:
        ok = 'output/share/alignments/contigs-to-ref/verkko_{release}_{commit}/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.copied.ok'
    run:
        import json
        import pathlib as pl
        import shutil as sh

        with open(input.version, 'r') as js_dump:
            verkko = json.load(js_dump)
            release = verkko['verkko_release']
            commit = verkko['verkko_commit']

        share_path = pl.Path(config['path_root_share_working']).resolve()
        verkko_subfolder = share_path / pl.Path(f'alignments/contigs-to-ref/verkko_{release}_{commit}/{wildcards.chrom}/{wildcards.hifi_type}')
        verkko_subfolder.mkdir(parents=True, exist_ok=True)

        check_file = ''

        for source in [input.bam, input.bai, input.paf]:
            source_path = pl.Path(source)
            dest_path = verkko_subfolder / source_path.name
            sh.copy(source_path, dest_path)
            check_file += f'{source_path}\t{dest_path}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK
