
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
