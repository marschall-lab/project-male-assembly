
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
        # motifs = expand(
        #     'references_derived/{motif}.fasta',
        #     motif=config['motif_search']
        # )
    output:
        ok = 'output/share/references.ok'
    run:
        import os
        import pathlib as pl
    
        share_path = pl.Path(config['path_root_share_references']).resolve()

        check_file = ''
        for ref_file in input:
            source = pl.Path(ref_file).resolve(strict=True)
            target = share_path / source.name
            rsync(source, target)
            check_file += f'{source}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_verkko_assemblies:
    """
    Collect additional files on-the-fly
    HPC outputs are not useful for downstream analyses (at the moment)
    """
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.info',
        linear = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
        index = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta.fai',
        hifi_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.hifi-coverage.csv',
        ont_cov = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.ont-coverage.csv',
    output:
        ok = 'output/share/assemblies/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    wildcard_constraints:
        mapq = 'na',
        chrom = 'wg'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'assemblies',
            wildcards.chrom
        )

        check_file = ''
        for source in [input.linear, input.index, input.hifi_cov, input.ont_cov]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_verkko_subset_assembly:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        linear = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.fasta',
        index = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.fasta.fai',
    output:
        ok = 'output/share/assemblies/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    wildcard_constraints:
        mapq = 'na',
        chrom = 'chrY'
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'assemblies',
            wildcards.chrom
        )

        check_file = ''
        for source in [input.linear, input.index]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_contig_to_ref_alignments:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.info',
        bam = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam',
        bai = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam.bai',
        paf = 'output/alignments/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz',
        bed = 'output/eval/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed',
        cov = 'output/eval/contigs-to-ref/10_renamed/{sample}.{hifi_type}.{ont_type}.na.wg_aln-to_{reference}.ref-cov.tsv',
    output:
        ok = 'output/share/alignments/contigs-to-ref/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.copied.ok'
    wildcard_constraints:
        mapq = 'na',
        chrom = 'wg'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/contigs-to-ref',
            wildcards.chrom
        )
        check_file = ''
        for source in [input.bam, input.bai, input.paf, input.bed, input.cov]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_subset_contig_to_ref_alignments:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        bam = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam',
        bai = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bam.bai',
        paf = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz',
    output:
        ok = 'output/share/alignments/contigs-to-ref/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.copied.ok'
    wildcard_constraints:
        mapq = 'na',
        chrom = 'chrY'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/contigs-to-ref',
            wildcards.chrom
        )
        check_file = ''
        for source in [input.bam, input.bai, input.paf]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


READ_ALN_SUBFOLDER = {
    'wg': 'alignments/reads-to-assm',
    'chrY': 'subset_wg/40_extract_rdaln'
}
rule copy_reads_to_assm_alignments:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        aln_files = lambda wildcards: expand(
            'output/{subfolder}/{{sample}}.{{other_reads}}_aln-to_{{hifi_type}}.{{ont_type}}.{{mapq}}.{{chrom}}.{ext}',
            subfolder=READ_ALN_SUBFOLDER[wildcards.chrom],
            ext=['bam', 'bam.bai', 'paf.gz']
        )
    output:
        ok = 'output/share/alignments/reads-to-assm/verkko_{major}_{minor}/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/reads-to-assm',
            wildcards.chrom
        )
        check_file = ''
        for source in input.aln_files:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_reads_to_ref_alignments:
    """
    So far, only requested for chrY
    """
    input:
        aln_files = lambda wildcards: expand(
            'output/subset_wg/60_subset_rdref/{{sample}}.{{other_reads}}_aln-to_{{reference}}.{{chrom}}.{ext}',
            ext=['bam', 'bam.bai', 'paf.gz']
        )
    output:
        ok = 'output/share/alignments/reads-to-ref/{sample}.{other_reads}_aln-to_{reference}.{chrom}.copied.ok'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl

        share_path = pl.Path(config['path_root_share_working']).resolve(strict=True)
        share_subfolder = share_path / pl.Path(f'alignments/reads-to-ref/{wildcards.chrom}')
        share_subfolder.mkdir(parents=True, exist_ok=True)

        check_file = ''
        for source in input.aln_files:
            source_path = pl.Path(source)
            target = share_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_motif_files:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        tsv = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
        bed = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm-hiq.bed',
        fasta = 'output/subset_wg/50_subset_motif/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hiq-seq.fasta',
    output:
        ok = 'output/share/motif_search/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.copied.ok'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'sequence_motifs',
            wildcards.chrom
        )
        check_file = ''
        for source in [input.bed, input.fasta, input.tsv]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_repeatmasker_files:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        fasta = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.rm-mask.fasta',
        table = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.matches.tsv',
        tar = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.rm-out.tar.gz',
    output:
        ok = 'output/share/repeatmasker/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'repeatmasker',
            wildcards.chrom
        )
        check_file = ''
        for source in [input.fasta, input.table, input.tar]:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_variant_calls:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        dv_vcf = 'output/variant_calls/10_filter_HIFIRW/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.dv-HET-SNV.vcf.gz',
        dv_tbi = 'output/variant_calls/10_filter_HIFIRW/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.dv-HET-SNV.vcf.gz.tbi',
        dv_stats = 'output/variant_calls/10_filter_HIFIRW/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.dv-HET-SNV.stats',
        pr_vcf = 'output/variant_calls/10_filter_ONTUL/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.pr-HET-SNV.vcf.gz',
        pr_tbi = 'output/variant_calls/10_filter_ONTUL/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.pr-HET-SNV.vcf.gz.tbi',
        pr_stats = 'output/variant_calls/10_filter_ONTUL/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.pr-HET-SNV.stats',
    output:
        ok = 'output/share/variant_calls/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'variant_calls',
            wildcards.chrom
        )
        check_file = ''
        source_files = [
            input.dv_vcf,
            input.dv_tbi,
            input.dv_stats,
            input.pr_vcf,
            input.pr_tbi,
            input.pr_stats
        ]
        for source in source_files:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


# DEPRECATED
#
# rule copy_chromosome_readsets:
#     input:
#         fasta = 'output/read_subsets/{chrom}/{sample}_{read_type}.{chrom}-reads.{mapq}.fasta.gz'
#     output:
#         ok = 'output/share/read_subsets/{sample}_{read_type}.{chrom}-reads.{mapq}.copied.ok'
#     run:
#         import os
#         import pathlib as pl
#         import shutil as sh

#         share_path = pl.Path(config['path_root_share_working']).resolve()
#         readsets_subfolder = share_path / pl.Path(f'readsets/{wildcards.chrom}/{wildcards.read_type}')
#         readsets_subfolder.mkdir(parents=True, exist_ok=True)

#         source = pl.Path(input.fasta)
#         destination = readsets_subfolder / source.name
#         if destination.is_file() and SKIP_SIZE:
#             source_size = os.stat(source).st_size
#             dest_size = os.stat(destination).st_size
#             if source_size != dest_size:
#                 sh.copy(source, destination)
#         else:
#             sh.copy(source, destination)

#         with open(output.ok, 'w') as dump:
#             _ = dump.write(f'{source}\t{destination}\n')
#     # END OF RUN BLOCK