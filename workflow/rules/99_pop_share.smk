
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
    priority: 100
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

        # EBI-UPLOAD
        # add whole-genome assemblies for EBI upload
        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'assemblies',
            "whole_genome",
            purpose="deposit"
        )
        for source in [input.linear, input.index]:
            source_path = pl.Path(source)
            sample = source_path.name.split(".")[0]
            if sample in QC_SAMPLES:
                target = verkko_subfolder.joinpath("qc_assembly", source_path.name)
            else:
                target = verkko_subfolder / source_path.name
            target.parent.mkdir(parents=True, exist_ok=True)
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
        chrom = '(chrY|chrX)'
    priority: 100
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

        # EBI-UPLOAD
        # add whole-genome assemblies for EBI upload
        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            "assemblies",
            wildcards.chrom,
            purpose="deposit"
        )
        for source in [input.linear, input.index]:
            source_path = pl.Path(source)
            sample = source_path.name.split(".")[0]
            if sample in QC_SAMPLES:
                target = verkko_subfolder.joinpath("qc_assembly", source_path.name)
            else:
                target = verkko_subfolder / source_path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_assembly_error_table:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        errors = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.errors.tsv',
    output:
        ok = 'output/share/assembly_errors/verkko_{major}_{minor}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    wildcard_constraints:
        mapq = 'na',
        chrom = 'chrY'
    priority: 100
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'assembly_errors',
            wildcards.chrom
        )

        check_file = ''
        for source in [input.errors]:
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
        chrom = '(chrY|chrX)'
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


rule copy_reference_motif_files:
    input:
        tsv = 'output/motif_search/10_norm/20_refseq/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
        bed = 'output/motif_search/10_norm/20_refseq/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm-hiq.bed',
        fasta = 'output/motif_search/10_norm/20_refseq/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hiq-seq.fasta',
    output:
        ok = 'output/share/reference_motifs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.copied.ok'
    wildcard_constraints:
        sample = 'T2T'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    run:
        import pathlib as pl

        subfolder = pl.Path(config['path_root_share_references']).resolve(strict=True)

        check_file = ''
        for source in [input.bed, input.fasta, input.tsv]:
            source_path = pl.Path(source)
            target = subfolder / source_path.name
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


rule copy_reference_repeatmasker_files:
    input:
        fasta = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.rm-mask.fasta',
        table = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.matches.tsv',
        tar = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.rm-out.tar.gz',
    output:
        ok = 'output/share/reference_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    wildcard_constraints:
        sample = 'T2T'
    resources:
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    run:
        import pathlib as pl

        subfolder = pl.Path(config['path_root_share_references']).resolve(strict=True)

        check_file = ''
        for source in [input.fasta, input.table, input.tar]:
            source_path = pl.Path(source)
            target = subfolder / source_path.name
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


rule copy_seq_class_alignments:
    input:
        version = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg.verkko.info',
        paf = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.paf.gz',
        pri = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.primary.bed',
        sec = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.secondary.bed',
        mrg = 'output/alignments/seqclasses-to-assm/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.merged.bed'
    output:
        ok = 'output/share/alignments/seqclasses-to-assm/verkko_{major}_{minor}/{seq_classes}_aln-to_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.copied.ok'
    run:
        import pathlib as pl

        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version,
            'alignments/seqclasses-to-assm',
            wildcards.chrom
        )
        check_file = ''
        source_files = [
            input.paf,
            input.pri,
            input.sec,
            input.mrg
        ]
        for source in source_files:
            source_path = pl.Path(source)
            target = verkko_subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_seqclass_annotations:
    input:
        version = expand(
            'output/hybrid/verkko/{sample}.HIFIRW.ONTUL.na.wg.verkko.info',
            sample=COMPLETE_SAMPLES
        ),
        t2t_bed = 'references_derived/T2T.chrY-seq-classes.bed',
        t2t_tsv = 'references_derived/T2T.chrY-seq-classes.tsv',
        t2t_fasta = 'references_derived/T2T.chrY-seq-classes.fasta',
        bed_generic = expand(
            'references_derived/seqclasses/{sample}.HIFIRW.ONTUL.na.chrY.generic-seqcls.bed',
            sample=COMPLETE_SAMPLES
        ),
        bed_specific = expand(
            'references_derived/seqclasses/{sample}.HIFIRW.ONTUL.na.chrY.specific-seqcls.bed',
            sample=COMPLETE_SAMPLES
        )
    output:
        ok = 'output/share/seqclasses/verkko_{major}_{minor}/ALL-SAMPLES_T2T-Y.copied.ok'
    run:
        import pathlib as pl
        import itertools as itt

        # EBI-UPLOAD only
        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version[0],
            "seq_classes",
            "all",
            purpose="deposit"
        )

        source_files = [
            input.t2t_bed, input.t2t_tsv, input.t2t_fasta,
        ]

        for source in itt.chain(source_files, input.bed_generic, input.bed_specific):
            source_path = pl.Path(source)
            sample = source_path.name.split(".")[0]
            if sample in QC_SAMPLES:
                target = verkko_subfolder.joinpath("qc_assembly", source_path.name)
            else:
                target = verkko_subfolder / source_path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_reference_motif_seq_files:
    input:
        motifs = expand(
            "references_derived/{motif}.fasta",
            motif=[
                "DYZ18_Yq", "DYZ19_Yq", "DYZ1_Yq",
                "DYZ2_Yq", "DYZ3-prim_Ycentro", "DYZ3-sec_Ycentro",
                "TSPY", "Yqhet_2k7bp", "Yqhet_3k1bp"
            ]
        )
    output:
        ok = "output/share/seq_motifs/ALL-MOTIFS.copied.ok"
    run:
        import pathlib as pl

        # EBI-UPLOAD only
        subfolder = pl.Path(config['path_root_deposit_ebi']).resolve(strict=True)
        subfolder = subfolder.joinpath("ref_seq_motifs")
        subfolder.mkdir(parents=True, exist_ok=True)

        check_file = ''
        for source in input.motifs:
            source_path = pl.Path(source)
            target = subfolder / source_path.name
            rsync(source_path, target)
            check_file += f'{source_path}\t{target}\n'

        with open(output.ok, 'w') as dump:
            _ = dump.write(check_file)
    # END OF RUN BLOCK


rule copy_flagged_regions:
    input:
        regions = expand(
            "output/eval/flagged_regions/sample_bed/{sample}.HIFIRW.ONTUL.na.chrY.flagged-all.bed",
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
        ),
        clusters = expand(
            "output/eval/flagged_regions/sample_bed/{sample}.HIFIRW.ONTUL.na.chrY.mixed-clusters.bed",
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
        )
    output:
        ok = 'output/share/flagged_regions/verkko_{major}_{minor}/ALL-SAMPLES_all-qc-tools.copied.ok',
    run:
        import pathlib as pl
        import itertools as itt

        # EBI-UPLOAD only
        verkko_subfolder, assembly_id = determine_verkko_subfolder(
            input.version[0],
            "flagged_regions",
            "all",
            purpose="deposit"
        )

        check_file = ''
        for source in itt.chain(input.regions, input.clusters):
            source_path = pl.Path(source)
            sample = source_path.name.split(".")[0]
            if sample in QC_SAMPLES:
                target = verkko_subfolder.joinpath("qc_assembly", source_path.name)
            else:
                target = verkko_subfolder / source_path.name
            target.parent.mkdir(parents=True, exist_ok=True)
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