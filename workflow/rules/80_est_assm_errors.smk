
########################################################################################
# Update 2022-05-25
# VerityMap (v2.1.1-alpha) can only be set up inside a Singularity container,
# because of an include error that appears when pip install is run
# in an (otherwise identical) Conda environment (see gh#18). However, the
# containerized version fails when processig wg assemblies w/o producing an
# informative error message, which makes debugging impossible. Hence,
# the ugly setup strategy below: build the environment in the container,
# but run setup outside of it to support code inspection (changes) on-the-fly.
# This creates the veritymap script in a .local cache folder, and other abominations...
# obvious TODO: replace that when errors are fixed
########################################################################################

localrules: clone_veritymap_repo, \
            build_veritymap, \
            clean_veritymap_build, \
            merge_error_annotations, \
            merge_all_error_stats, \
            aggregate_errors_per_seqclass, \
            merge_agg_seqclass_errors


rule clone_veritymap_repo:
    """
    2022-07-04
    Checkout target is bugfix branch
    """
    output:
        'repos/veritymap.git.ok'
    envmodules:
        'git/2.24.0'
    params:
        veritymap_repo = config['veritymap_repo']
    shell:
        'mkdir -p repos/ && cd repos/'
            ' && '
        'git clone {params.veritymap_repo}'
            ' && '
        'cd VerityMap'
            ' && '
        'git checkout 8d241f4'
            ' && '
        'cd ../../ && touch {output}'


rule clean_veritymap_build:
    input:
        'repos/veritymap.git.ok'
    output:
        'repos/veritymap.clean.ok'
    shell:
        'rm -rfd .cache/pip'
            ' && '
        'rm -rfd .local/lib && rm -rfd .local/bin'
            ' && '
        'rm -rfd repos/VerityMap/.cache'
            ' && '
        'rm -rfd repos/VerityMap/VerityMap.egg-info'
            ' && '
        'rm -rfd repos/VerityMap/build'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/build'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/__pycache__'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/py_src/__pycache__'
            ' && '
        'touch {output}'


rule build_veritymap:
    input:
        'repos/veritymap.git.ok',
        'repos/veritymap.clean.ok'
    output:
        'repos/veritymap.build.ok'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap_env"]}'
    shell:
        'pip install --no-warn-script-location --no-deps repos/VerityMap/'
            ' && '
        'touch {output}'


VERITYMAP_MEM_FACTOR = {
    'ONTUL': 4
}


rule run_chrom_veritymap:
    """
    NB: important to execute the main.py at the right
    location due to the setup horror explained above.
    Since main is messing with sys.path, VerityMap
    can be properly executed irrespective of the install/setup
    into "." locations by pip.
    """
    input:
        build_ok = 'repos/veritymap.build.ok',
        reads = 'output/subset_wg/45_extract_reads/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.reads.fasta.gz',
        assembly = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.fasta',
    output:
        chk = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.chk',
        bed = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}/{sample}_kmers_dist_diff.bed'
    log:
        'log/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.log'
    benchmark:
        'rsrc/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.rsrc'
    wildcard_constraints:
        chrom = '(chrY|chrX)'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap_env"]}'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt * VERITYMAP_MEM_FACTOR.get(wildcards.other_reads, 1),
        walltime = lambda wildcards, attempt: f'{11*attempt}:59:00'
    params:
        preset = lambda wildcards: {'HIFIRW': 'hifi', 'ONTUL': 'ont'}[wildcards.other_reads],
        outdir = lambda wildcards, output: output.chk.rsplit('.', 2)[0]
    shell:
        'python repos/VerityMap/veritymap/main.py --no-reuse --reads {input.reads} -t {threads} '
        '-d {params.preset} -l {wildcards.sample} -o {params.outdir} {input.assembly} &> {log}'
            ' && '
        'touch {output.chk}'


rule normalize_veritymap_bed_file:
    """
    Normalize VerityMap BED file for simple
    merge with SNV/HET errors from variant
    callers - just format change...

    2022-07-08
    Feedback from dev indicates that records that
    all give the same misassembly length but different
    coordinates for the rare k-mers are - more likely
    than not - all describing the same problem/mis-
    assembly event, and VerityMap is just internally
    failing at recognizing this. Not waiting for a
    fix, but merging all such records into one keeping
    lowest start and highest end coordinate.
    """
    input:
        bed = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}/{sample}_kmers_dist_diff.bed'
    output:
        tsv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm-errors.tsv'
    run:
        import pandas as pd

        keep_cols = ['chrom', 'start', 'end', 'est_size']
        df = pd.read_csv(
            input.bed, sep='\t', header=None,
            names=keep_cols + ['num_reads'],
            usecols=keep_cols
        )
        df.sort_values(['chrom', 'start'], ascending=True, inplace=True)

        # 2022-07-08 (see docstring above)
        merged_events = []
        for (ctg, misassm_len), records in df.groupby(['chrom', 'est_size']):
            if records.shape[0] == 1:
                # unique event
                merged_events.append(records.copy())
                continue
            
            first_iv_start, first_iv_end = records.loc[records.index[0], ['start', 'end']].values
            last_iv_start, last_iv_end = records.loc[records.index[-1], ['start', 'end']].values

            # heuristic 1: if all the same event,
            # the estimated size should span all intervals
            is_spanning = (last_iv_end - first_iv_start) <= abs(misassm_len)

            # heuristic 2: the event should at least
            # reach into first and last interval
            is_reaching = (first_iv_end - 1 + misassm_len) >= last_iv_start

            # heuristic 3: more like a bug fix heuristic;
            # in some cases, VerityMap seems to mess up
            # coordinates and reports zero-length
            # region(s) instead of a single one with L > 0
            any_zero_length = (records.shape[0] == 2) & ((records['end'] - records['start']) == 0).any()

            # heuristic 4: if the event is small enough
            # to actually be contained in all reported
            # regions, then assume the events are indeed
            # disjoint (could be a problem in repeat
            # resolution)
            is_disjoint = (records['end'] - records['start'] > abs(misassm_len)).all()

            if is_spanning or is_reaching or any_zero_length:
                new_record = pd.DataFrame(
                    [[ctg, first_iv_start, last_iv_end, misassm_len]],
                    columns=keep_cols
                )
                merged_events.append(new_record)
            elif is_disjoint:
                merged_events.append(records.copy())
            else:
                # What is left at this point: records that are not clearly
                # disjoint but also not spanning or reaching - be conservative
                # and create separate "reaching" records. In essence, this will
                # probably overestimate the actual error rate.
                cover_reach = abs(misassm_len)
                covered_indices = set()
                for row in records.itertuples(index=False):
                    start = row.start
                    end = start + cover_reach
                    select_reach = (records['start'] < end) & (start < records['end'])
                    sub = records.loc[select_reach, :]
                    if all(i in covered_indices for i in sub.index):
                        continue
                    new_record = pd.DataFrame(
                        [[ctg, sub['start'].min(), sub['end'].max(), misassm_len]],
                        columns=keep_cols
                    )
                    merged_events.append(new_record)
                    covered_indices = covered_indices.union(set(sub.index.values))
                    
                #raise ValueError(f'Cannot process presumably identical events: {records}')

        merged_events = pd.concat(merged_events, axis=0, ignore_index=False)    

        merged_events['sample'] = wildcards.sample
        merged_events['reads'] = wildcards.other_reads
        merged_events['source'] = 'VerityMap'
        merged_events.sort_values(['chrom', 'start'], ascending=True, inplace=True)

        assert pd.notnull(merged_events).all(axis=1).all()

        merged_events.to_csv(output.tsv, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule merge_error_annotations:
    input:
        snv_dv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-errors.tsv',
        snv_pr = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-errors.tsv',
        hifi_vm = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-errors.tsv',
        quast = 'output/eval/assm_stats/SAMPLES.{hifi_type}.{ont_type}.na.{chrom}.quast-report.tsv'
    output:
        tsv = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.errors.tsv',
        stats = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.error-stats.tsv',
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    run:
        import pandas as pd

        assm_stats = pd.read_csv(input.quast, sep='\t', header=0)
        assm_length = assm_stats.loc[assm_stats['sample'] == wildcards.sample, 'assembly_length_bp'].values[0]

        all_errors = []
        error_stats = dict()
        for err_file in [input.snv_dv, input.snv_pr, input.hifi_vm]:
            df = pd.read_csv(err_file, sep='\t', header=0)
            all_errors.append(df)
        all_errors = pd.concat(all_errors, axis=0, ignore_index=False)
        all_errors.sort_values(['chrom', 'start', 'end'], ascending=True, inplace=True)

        error_stats['assembly_length_bp'] = assm_length
        error_stats['num_errors'] = all_errors.shape[0]
        error_stats['num_SNV_errors'] = all_errors.loc[all_errors['est_size'] == 1, :].shape[0]
        error_stats['error_bp'] = all_errors['est_size'].abs().sum()
        error_stats['error_bp_per_kbp'] = round(error_stats['error_bp'] / assm_length * 1000, 5)
        error_stats['error_SNV_per_kbp'] = round(error_stats['num_SNV_errors'] / assm_length * 1000, 5)
        error_stats['sample'] = wildcards.sample

        all_errors.to_csv(output.tsv, sep='\t', header=True, index=False)

        error_stats = pd.DataFrame.from_records([error_stats])
        error_stats.to_csv(output.stats, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule merge_all_error_stats:
    input:
        stats = expand(
            'output/eval/merged_errors/{sample}.{{hifi_type}}.{{ont_type}}.na.{{chrom}}.error-stats.tsv',
            sample=COMPLETE_SAMPLES
        )
    output:
        tsv = 'output/eval/merged_errors/SAMPLES.{hifi_type}.{ont_type}.na.{chrom}.error-stats.tsv',
    run:
        import pandas as pd

        all_stats = []
        for tsv_file in input.stats:
            df = pd.read_csv(tsv_file, sep='\t', header=0)
            all_stats.append(df)

        all_stats = pd.concat(all_stats, axis=0, ignore_index=False)
        all_stats['qc_only_assembly'] = 0
        qc_only = all_stats['sample'].isin(QC_SAMPLES)
        all_stats.loc[qc_only, 'qc_only_assembly'] = 1

        all_stats.sort_values(['qc_only_assembly', 'sample'], ascending=True, inplace=True)

        with open(output.tsv, 'w') as table:
            _ = table.write('## 80_est_assm_errors::merge_all_error_stats\n')
            _ = table.write('## Suppl. table listing error summary statistics per assembly\n')

        all_stats.to_csv(output.tsv, sep='\t', mode='a', header=True, index=False)
    # END OF RUN BLOCK


rule intersect_errors_with_seqclasses:
    """
    2022-07-27
    After updating the preprocessing for the sequence class annotation,
    the input BED has been simplified to just contain plain seq. class
    names (same as for T2T-Y reference)
    """
    input:
        errors = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.errors.tsv',
        seqclasses = 'references_derived/seqclasses/{sample}.{hifi_type}.{ont_type}.na.{chrom}.generic-seqcls.bed',
    output:
        table = 'output/eval/error_clusters/10_intersect/{sample}.{hifi_type}.{ont_type}.na.{chrom}.isect.tsv',
        tmp = temp('output/eval/error_clusters/{sample}.{hifi_type}.{ont_type}.na.{chrom}.errors.bed')
    wildcard_constraints:
        chrom = 'chrY'
    conda:
        '../envs/biotools.yaml'
    shell:
        'tail -n +2 {input.errors} > {output.tmp}'
            ' && '
        'bedtools intersect -wao -a {input.seqclasses} -b {output.tmp} > {output.table}'


rule aggregate_errors_per_seqclass:
    """
    2022-07-27
    See update above; consequently, normalization of region names is
    no longer needed (removed internal funtion "norm_region_name")
    """
    input:
        table = 'output/eval/error_clusters/10_intersect/{sample}.{hifi_type}.{ont_type}.na.{chrom}.isect.tsv'
    output:
        table = 'output/eval/error_clusters/20_aggregate/{sample}.{hifi_type}.{ont_type}.na.{chrom}.agg-seqclass-errors.tsv'
    run:
        import pandas as pd
        import re as re

        columns = [
            'contig', 'start', 'end', 'region_type',
            'contig2', 'err_start', 'err_end', 'err_size',
            'sample', 'err_source', 'err_reads', 'overlap'
        ]

        df = pd.read_csv(input.table, sep='\t', header=None, names=columns)
        df['err_size'].replace('.', '0', inplace=True)
        df['err_size'] = df['err_size'].astype(int)
        df['err_size'] = df['err_size'].abs()  # DEL are given as neg.
        df['sample'] = wildcards.sample  # no-hits in the intersection have "." as sample

        records = []
        for (ctg, region), errors in df.groupby(['contig', 'region_type']):
            sample = errors.at[errors.index[0], 'sample']
            assert sample in ctg
            total_err_bp = errors['err_size'].sum()
            total_err_num = errors.shape[0]
            if total_err_bp == 0:
                total_err_num = 0
            region_size = errors.at[errors.index[0], 'end'] - errors.at[errors.index[0], 'start']
            assert region_size > 0
            records.append((sample, ctg, region, region_size, total_err_bp, total_err_num))

        out_columns = ['sample', 'contig', 'region_type', 'region_size', 'errors_bp', 'errors_num']
        split_stats = pd.DataFrame.from_records(records, columns=out_columns)

        agg_stats = split_stats.groupby(['sample', 'region_type']).agg(
            contigs_num=('contig', 'nunique'),
            errors_bp=('errors_bp', sum),
            errors_num=('errors_num', sum),
            region_size=('region_size', sum)
        )
        agg_stats.reset_index(drop=False, inplace=True)
        agg_stats.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule merge_agg_seqclass_errors:
    input:
        tables = expand(
            'output/eval/error_clusters/20_aggregate/{sample}.{{hifi_type}}.{{ont_type}}.{{mapq}}.{{chrom}}.agg-seqclass-errors.tsv',
            sample=[s for s in COMPLETE_SAMPLES if s != 'NA24385']
        )
    output:
        table = 'output/eval/error_clusters/SAMPLES.{hifi_type}.{ont_type}.{mapq}.{chrom}.mrg-seqclass-errors.tsv'
    run:
        import pandas as pd
        merged = []
        for table in input.tables:
            df = pd.read_csv(table, sep='\t', header=0)
            merged.append(df)
        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged['qc_only_assembly'] = 0
        qc_only = merged['sample'].isin(QC_SAMPLES)
        merged.loc[qc_only, 'qc_only_assembly'] = 1
        merged.sort_values(['qc_only_assembly', 'sample', 'region_type'], ascending=True, inplace=True)
        merged.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


###################################################
## Below: new rules introduced for revised version
###################################################


rule filter_chromosome_alignments:
    """This rule filters the read-to-assembly
    alignments (only filter: discard unmapped)
    to remove secondary and supplementary
    alignments following the README of
    NucFreq: https://github.com/mrvollger/NucFreq
        > we filtered the alignments to remove secondary
        > and partial alignment using SAMtools flag 2308
    """
    input:
        bam = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.bam.bai',
    output:
        bam = 'output/subset_wg/40_extract_rdaln/nucfreq/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.filtered.bam',
    conda:
        "../envs/biotools.yaml"
    threads: config["num_cpu_low"]
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        sam_flag = 2308
    shell:
        "samtools view -F {params.sam_flag} --threads {threads} "
        "--bam -o {output.bam} {input.bam}"


localrules: filter_contigs_for_nucfreq
rule filter_contigs_for_nucfreq:
    input:
        bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
    output:
        bed = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.ctg-500kbp.bed',
    wildcard_constraints:
        chrom = '(chrY|chrX)'
    run:
        import pandas as pd
        names = ["contig", "start", "end"]
        contigs = pd.read_csv(input.bed, sep="\t", header=None, names=names)
        contigs = contigs.loc[contigs["end"] >= 500000, :].copy()
        contigs.sort_values("contig", ascending=True, inplace=True)
        contigs.to_csv(output.bed, sep="\t", header=False, index=False)
    # END OF RUN BLOCK


rule run_nucfreq:
    input:
        bam = 'output/subset_wg/40_extract_rdaln/nucfreq/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.filtered.bam',
        bai = 'output/subset_wg/40_extract_rdaln/nucfreq/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.chrY.filtered.bam.bai',
        bed = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.ctg-500kbp.bed',
    output:
        bed = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.nucfreq.bed',
        png = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.nucfreq.png'
    log:
        'log/output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.nucfreq.log'
    singularity:
        f"{config['container_store']}/{config['container']['nucfreq']}"
    threads: config["num_cpu_medium"]
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt}:59:00',
    shell:
        "NucPlot.py --obed {output.bed} --threads {threads} --bed {input.bed} "
        "{input.bam} {output.png} &> {log}"


rule zip_nucfreq_plots:
    """ NucFreq errors out
    for HG00512
    """
    input:
        png = expand(
            'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{{chrom}}.{other_reads}.nucfreq.png',
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            other_reads=["HIFIRW"]
        )
    output:
        zipped = 'output/eval/assm_errors/nucfreq/ALL-SAMPLES.{chrom}.nucfreq-plots.zip',
    shell:
        "zip -9 -D -j {output.zipped} {input.png}"



localrules: detect_nucfreq_het_positions
rule detect_nucfreq_het_positions:
    """Following descriptions in NucFreq README:

        > identify regions where the second most
        > common base was present in at least 10%
        > of reads in at least 5 positions within
        > a 500 bp region.
    """
    input:
        bed = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.nucfreq.bed',
    output:
        het = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.het-regions.tsv',
    run:
        import pandas as pd
        names = ["contig", "start", "end", "first", "second"]
        df = pd.read_csv(input.bed, sep="\t", header=None, comment="#", names=names)

        # het_ratio will be aggregated with median later,
        # hence naming here like this to avoid renaming
        df["median_het_ratio"] = (df["second"]/(df["first"]+df["second"]) * 100).round(1)
        # (1) present in at least 10%
        df = df.loc[df["median_het_ratio"] >= 10, :].copy()
        # (2) within a 500 bp region
        # NB: use interval end here to
        # use max HET position in aggregate
        # downstream (that is typically smaller)
        df["iv_end"] = df["start"] + 500
        # need the following indicator to track
        # how many rows were merged per interval
        df["num_hets"] = 1

        hets_per_seq = []
        for contig, hets in df.groupby("contig"):
            hets["iv_idx"] = (hets["start"]>hets["iv_end"].shift().cummax()).cumsum()
            flag_regions = hets.groupby("iv_idx").agg(
                {
                    "start":"min", "end": "max",
                    "num_hets": "sum",
                    "median_het_ratio": "median"  # here: hence the name median_het_ratio
                }
            )
            # (3) at least 5 het positions
            flag_regions = flag_regions.loc[flag_regions["num_hets"] > 4, :].copy()
            flag_regions["contig"] = contig
            hets_per_seq.append(flag_regions)
            
        hets = pd.concat(hets_per_seq, axis=0, ignore_index=False)
        hets.sort_values(["contig", "start"], ascending=True, inplace=True)
        hets = hets[["contig", "start", "end", "median_het_ratio", "num_hets"]]
        hets.reset_index(drop=True, inplace=True)
        hets.to_csv(output.het, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: compute_nucfreq_sample_stats
rule compute_nucfreq_sample_stats:
    input:
        seqclasses = 'references_derived/seqclasses/{sample}.{hifi_type}.{ont_type}.na.{chrom}.generic-seqcls.bed',
        chrom_bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
        ctg_bed = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.ctg-500kbp.bed',
        het_regions = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.het-regions.tsv',
    output:
        stats = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.stats.tsv'
    run:
        import pandas as pd

        def compute_genome_size(fp):
            df = pd.read_csv(fp, sep="\t", header=None, names=["contig", "start", "end"])
            df["length"] = df["end"] - df["start"]
            genome_size = int(df["length"].sum())
            num_contigs = int(df["contig"].nunique())
            return num_contigs, genome_size

        def get_het_seqclasses(fp):
            seqcls = pd.read_csv(
                fp, sep="\t", header=None,
                names=["contig", "start", "end", "seqclass"]
            )
            # HET3_Yq
            seqcls["is_YqHET"] = seqcls["seqclass"].apply(lambda x: 1 if "HET3" in x else 0)
            assert seqcls["is_YqHET"].sum() > 0
            # HET1_centro
            seqcls["is_CEN"] = seqcls["seqclass"].apply(lambda x: 1 if "HET1" in x else 0)
            assert seqcls["is_CEN"].sum() > 0
            get_yqhet = seqcls["is_YqHET"] == 1
            get_cen = seqcls["is_CEN"] == 1
            seqcls = seqcls.loc[get_yqhet | get_cen, :].copy()
            return seqcls

        def find_het_region_ovl_het(het_seqclasses, het_regions):
            region_in_yqhet = 0
            region_in_cen = 0
            for ctg, regions in het_regions.groupby("contig"):
                classes = het_seqclasses.loc[het_seqclasses["contig"] == ctg, :]
                for region in regions.itertuples(index=True):
                    if region.end < classes["start"].min() or region.start > classes["end"].max():
                        continue
                    get_start = region.start < classes["end"]
                    get_end = region.end > classes["start"]
                    get_ovl = get_start & get_end
                    overlapping = classes.loc[get_ovl, :]
                    if overlapping.empty:
                        continue
                    if overlapping["is_YqHET"].sum() > 0:
                        region_in_yqhet += 1
                    if overlapping["is_CEN"].sum() > 0:
                        region_in_yqhet += 1

            num_regions = het_regions.shape[0]
            region_stats = {
                "num_HET_regions": num_regions,
                "num_ovl_YqHET": region_in_yqhet,
                "num_ovl_CEN": region_in_cen,
                "num_ovl_YqHET_or_CEN": region_in_yqhet + region_in_cen,
                "pct_ovl_YqHET": round(region_in_yqhet / num_regions * 100, 1),
                "pct_ovl_CEN": round(region_in_cen / num_regions * 100, 1),
                "pct_ovl_YqHET_or_CEN": round((region_in_cen + region_in_yqhet) / num_regions * 100, 1),
            }
            return region_stats

        contigs_chrom, size_chrom = compute_genome_size(input.chrom_bed)
        contigs_500, size_500 = compute_genome_size(input.ctg_bed)
        size_pct = round(size_500 / size_chrom * 100, 2)

        sample_stats = {
            "sample": wildcards.sample,
            f"contigs_num_{wildcards.chrom}": contigs_chrom,
            f"size_bp_{wildcards.chrom}": size_chrom,
            f"contigs_num_{wildcards.chrom}_geq500kbp": contigs_500,
            f"size_bp_{wildcards.chrom}_geq500kbp": size_500,
            f"processed_{wildcards.chrom}_geq500kbp": size_pct
        }
        df = pd.read_csv(input.het_regions, sep="\t", header=0)
        region_stats = find_het_region_ovl_het(
            get_het_seqclasses(input.seqclasses),
            df
        )
        sample_stats.update(region_stats)

        df["het_region_length"] = df["end"] - df["start"]
        agg = df.agg(
            contigs_num=("contig", "nunique"),
            total_flagged_length_bp=("het_region_length", sum),
            total_hets_num=("num_hets", sum),
            min_agg_het_ratio=("median_het_ratio", min),
            max_agg_het_ratio=("median_het_ratio", max),
            median_agg_het_ratio=("median_het_ratio", "median"),
        )

        for label in agg.index:
            selector = pd.isnull(agg.loc[label, :].values)
            value = agg.loc[label, ~selector]
            if label.endswith("ratio"):
                value = float(value)
            else:
                value = int(value)
            sample_stats[label] = value
        sample_stats = pd.DataFrame.from_records([sample_stats])
        sample_stats.to_csv(output.stats, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: merge_all_nucfreq_stats
rule merge_all_nucfreq_stats:
    input:
        tables = expand(
            'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{{chrom}}.{other_reads}.stats.tsv',
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            other_reads=["HIFIRW"]
        )
    output:
        table = 'output/eval/assm_errors/nucfreq/ALL-SAMPLES.{chrom}.nucfreq-stats.tsv',
    run:
        import pandas as pd

        merged = []
        for table in input.tables:
            df = pd.read_csv(table, sep="\t", header=0)
            df["is_qc_sample"] = 0
            if df["sample"].values[0] in QC_SAMPLES:
                df["is_qc_sample"] = 1
            merged.append(df)

        merged = pd.concat(merged, axis=0, ignore_index=False)
        size_column = f"size_bp_{wildcards.chrom}_geq500kbp"
        merged["het_per_kbp"] = (merged["total_hets_num"] / (merged[size_column] / 1e3)).round(4)
        merged.sort_values(["is_qc_sample", "sample"], inplace=True)
        merged.to_csv(output.table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


############################################################
# reprocess all flagged regions to check for mutual support
# in a more consistent manner

rule convert_to_flanked_genomic_region:
    input:
        snv_dv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-errors.tsv',
        snv_pr = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-errors.tsv',
        vm_reg = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-errors.tsv',
        nf_reg = 'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.het-regions.tsv',        
        labels = 'references_derived/veritymap_expert-label.tsv',
    output:
        table_dv = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-genreg.tsv',
        table_pr = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-genreg.tsv',
        table_vm = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-genreg.tsv',
        table_nf = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.nf-genreg.tsv',
        bed_dv = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-genreg.bed',
        bed_pr = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-genreg.bed',
        bed_vm = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-genreg.bed',
        bed_nf = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.nf-genreg.bed',
    conda:
        '../envs/pyscript.yaml'
    params:
        script = find_script_path("norm_flagged_regions.py")
    shell:
        "{params.script} --table {input.snv_dv} --output {output.table_dv}"
            " && "
        "{params.script} --table {input.snv_pr} --output {output.table_pr}"
            " && "
        "{params.script} --table {input.vm_reg} --output {output.table_vm} --curated {input.labels}"
            " && "
        "{params.script} --table {input.nf_reg} --output {output.table_nf}"


rule read_depth_in_flagged_regions:
    """
    NB: this is an adapted c&p of the rule:
    63_eval_read_align::dump_read_to_assm_coverage

    NB: samtools depth skips reads ...

        > By default, reads that have any of the flags
        > UNMAP, SECONDARY, QCFAIL, or DUP set are skipped.
        (see man pages)

    """
    input:
        bed = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.{tool}-genreg.bed',
        bam = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam',
        bai = 'output/alignments/reads-to-assm/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.wg.bam.bai'
    output:
        depth = 'output/eval/flagged_regions/read_depth/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.minmq{minmapq}.{tool}-genreg.depth.tsv.gz'
    wildcard_constraints:
        other_reads = "(HIFIRW|ONTUL)",
        tool = "(vm|nf)"  # VerityMap or NucFreq
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt}:59:00'
    shell:
        'samtools depth -a --min-MQ {wildcards.minmapq} -l 5000 '
        '-b {input.bed} {input.bam} | pigz -p 2 > {output.depth}'


rule read_depth_in_chrom_y:
    """
    NB: this is an adapted c&p of the rule:
    63_eval_read_align::dump_read_to_assm_coverage

    NB: samtools depth skips reads ...

        > By default, reads that have any of the flags
        > UNMAP, SECONDARY, QCFAIL, or DUP set are skipped.
        (see man pages)

    """
    input:
        bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
        bam = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.bam',
        bai = 'output/subset_wg/40_extract_rdaln/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.bam.bai',
    output:
        cov = 'output/eval/flagged_regions/chrom_cov/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.minmq{minmapq}.cov.tsv.gz'
    wildcard_constraints:
        other_reads = "(HIFIRW|ONTUL)",
        chrom = "chrY"
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt}:59:00'
    shell:
        'samtools depth -a --min-MQ {wildcards.minmapq} -l 5000 '
        '-b {input.bed} {input.bam} | cut -f 3 | pigz -p 2 > {output.cov}'



rule cluster_flagged_regions:
    input:
        bed_dv = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-genreg.bed',
        bed_pr = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-genreg.bed',
        bed_vm = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-genreg.bed',
        bed_nf = 'output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.nf-genreg.bed',
    output:
        concat = temp('output/eval/flagged_regions/clustered/{sample}.{hifi_type}.{ont_type}.na.{chrom}.concat.bed'),
        merged = 'output/eval/flagged_regions/clustered/{sample}.{hifi_type}.{ont_type}.na.{chrom}.cluster.tsv',
    conda:
        '../envs/biotools.yaml'
    shell:
        "sort -V -k1,1 -k2n,3n {input} > {output.concat}"
            " && "
        "bedtools merge -c 4 -o collapse -i {output.concat} > {output.merged}"


rule identify_mixed_support_clusters:
    input:
        merged = "output/eval/flagged_regions/clustered/{sample}.{hifi_type}.{ont_type}.na.{chrom}.cluster.tsv",
        depths = expand(
            "output/eval/flagged_regions/read_depth/{{sample}}.{other_reads}_aln-to_{{hifi_type}}.{{ont_type}}.na.{{chrom}}.minmq{minmapq}.{tool}-genreg.depth.tsv.gz",
            other_reads=["HIFIRW", "ONTUL"],
            tool=["vm", "nf"],
            minmapq=[0, 10]
        ),
        coverages = expand(
            "output/eval/flagged_regions/chrom_cov/{{sample}}.{other_reads}_aln-to_{{hifi_type}}.{{ont_type}}.na.{{chrom}}.minmq{minmapq}.cov.tsv.gz",
            other_reads=["HIFIRW", "ONTUL"],
            minmapq=[0, 10]
        ),
        hifi_regions = expand(
            "output/eval/flagged_regions/flanked/{{sample}}.{{hifi_type}}.{{ont_type}}.na.{{chrom}}.HIFIRW.{tool}-genreg.tsv",
            tool=["vm", "nf", "dv"]
        ),
        ont_regions = "output/eval/flagged_regions/flanked/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-genreg.tsv"
    output:
        tsv = "output/eval/flagged_regions/annotated/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-clustered.tsv"
    conda:
        '../envs/pyscript.yaml'
    params:
        script = find_script_path("define_support_clusters.py")
    shell:
        "{params.script} --regions {input.hifi_regions} {input.ont_regions} "
            "--coverages {input.coverages} --depths {input.depths} "
            "--merged {input.merged} --output {output.tsv}"


rule dump_sample_beds_flagged_regions:
    input:
        tsv = "output/eval/flagged_regions/annotated/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-clustered.tsv",
        chrom_bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
    output:
        bed_regions = "output/eval/flagged_regions/sample_bed/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-all.bed",
        bed_clusters = "output/eval/flagged_regions/sample_bed/{sample}.{hifi_type}.{ont_type}.na.{chrom}.mixed-clusters.bed",
    run:
        import pandas as pd

        regions = pd.read_csv(input.tsv, header=0, sep="\t")
        regions = regions.loc[regions["region_type"] == "origin", :].copy()
        regions.sort_values(["contig", "start", "end"], inplace=True)
        regions["is_cluster_member"] = regions["cluster_type"].apply(lambda x: 0 if x == "singleton" else 1)
        bed_regions_columns = [
            "contig", "start", "end",
            "name", "software", "is_cluster_member",
            "cluster_id", "cluster_type"
        ]
        with open(output.bed_regions, "w") as dump:
            _ = dump.write("#")
            regions[bed_regions_columns].to_csv(dump, header=True, index=False, sep="\t")

        bed_cluster_columns = [
            "contig", "cluster_start", "cluster_end", "cluster_id",
            "num_nucfreq_regions", "num_veritymap_regions",
            "num_het_snv_hifi", "num_het_snv_ont", "snv_density_kbp"
        ]

        with open(output.bed_clusters, "w") as dump:
            _ = dump.write("#")
            regions.loc[regions["cluster_type"] == "mixed_regions", bed_cluster_columns].to_csv(
                dump, header=True, index=False, sep="\t"
            )
    # END OF RUN BLOCK


localrules: dump_sample_stats_flagged_regions
rule dump_sample_stats_flagged_regions:
    input:
        tsv = "output/eval/flagged_regions/annotated/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-clustered.tsv",
        chrom_bed = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.bed',
    output:
        sample_stats = "output/eval/flagged_regions/sample_stats/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-stats.tsv"
    run:
        import pandas as pd
        import collections as col

        def compute_genome_size(fp):
            df = pd.read_csv(fp, sep="\t", header=None, names=["contig", "start", "end"])
            df["length"] = df["end"] - df["start"]
            genome_size = int(df["length"].sum())
            num_contigs = int(df["contig"].nunique())
            return num_contigs, genome_size

        contigs, wg_size = compute_genome_size(input.chrom_bed)

        regions = pd.read_csv(input.tsv, header=0, sep="\t")
        regions = regions.loc[regions["region_type"] == "origin", :].copy()
        regions["length"] = regions["end"] - regions["start"]

        select_nucfreq = regions["software"] == "NucFreq"
        select_veritymap = regions["software"] == "VerityMap"
        select_pos = regions["software"].isin(["DeepVariant", "PEPPER"])
        select_reg = regions["software"].isin(["NucFreq", "VerityMap"])

        clusters = regions.loc[regions["cluster_type"] == "mixed_regions", :].copy()
        clusters.drop_duplicates("cluster_id", inplace=True, keep="first")

        stats = col.OrderedDict([
            ("sample", wildcards.sample),
            ("mixed_region_clusters_num", clusters.shape[0]),
            ("mixed_region_clusters_bp", int(clusters["cluster_span"].sum())),
            ("mixed_region_clusters_pct", 0),
            ("mixed_region_clusters_median_size", 0),
            ("mixed_region_clusters_mean_size", 0),
            ("flagged_regions_num", int(regions.loc[select_reg, :].shape[0])),
            ("flagged_regions_bp", int(regions.loc[select_reg, "length"].sum())),
            ("flagged_regions_pct", 0),
            ("het_snv_num", int(regions.loc[select_pos, "length"].sum())),
            ("het_snv_kbp_density", 0),
            ("flagged_nucfreq_num", int(regions.loc[select_nucfreq, :].shape[0])),
            ("flagged_nucfreq_bp", int(regions.loc[select_nucfreq, "length"].sum())),
            ("flagged_veritymap_num", int(regions.loc[select_veritymap, :].shape[0])),
            ("flagged_veritymap_bp", int(regions.loc[select_veritymap, "length"].sum()))
        ])

        if stats["flagged_regions_bp"] > 0:
            pct_flagged = round(stats["flagged_regions_bp"] / wg_size * 100, 2)
            stats["flagged_regions_pct"] = pct_flagged

        if stats["het_snv_num"] > 0:
            snv_density = round(stats["het_snv_num"] / (wg_size / 1000), 2)
            stats["het_snv_kbp_density"] = snv_density

        if stats["mixed_region_clusters_bp"] > 0:
            clustered_pct = round(stats["mixed_region_clusters_bp"] / wg_size * 100, 2)
            stats["mixed_region_clusters_pct"] = clustered_pct      

        if stats["mixed_region_clusters_num"] > 0:
            stats["mixed_region_clusters_median_size"] = int(clusters["cluster_span"].median())
            stats["mixed_region_clusters_mean_size"] = int(clusters["cluster_span"].mean())


        df = pd.DataFrame.from_records([stats])
        df.to_csv(output.sample_stats, header=True, index=False, sep="\t")
    # END OF RUN BLOCK


rule run_all_assm_errors:
    """
    Outlier: HG00512 is too fragmented
    and generates a plot output that is
    too large and cannot be created
    - remove for now
    """
    input:
        nucfreq = expand(
            'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.het-regions.tsv',
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            chrom=["chrY"],
            other_reads=["HIFIRW"]
        ),
        nucfreq_plots = expand(
            'output/eval/assm_errors/nucfreq/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.nucfreq.png',
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            chrom=["chrY"],
            other_reads=["HIFIRW"]
        ),
        zipped_plots = expand(
            'output/eval/assm_errors/nucfreq/ALL-SAMPLES.{chrom}.nucfreq-plots.zip',
            chrom=["chrY"]
        ),
        stats_table = expand(
            'output/eval/assm_errors/nucfreq/ALL-SAMPLES.{chrom}.nucfreq-stats.tsv',
            chrom=["chrY"]
        ),
        flagged_stats = expand(
            "output/eval/flagged_regions/sample_stats/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-stats.tsv",
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            chrom=["chrY"],
        ),
        flagged_bed = expand(
            "output/eval/flagged_regions/sample_bed/{sample}.{hifi_type}.{ont_type}.na.{chrom}.flagged-all.bed",
            sample=[s for s in COMPLETE_SAMPLES if s != "HG00512"],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            chrom=["chrY"],
        )
