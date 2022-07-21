
rule one_hot_encoding_projects_haplogroups:
    input:
        table = config['projects']
    output:
        table = 'output/quality_factors/project_haplogroup_feature.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep='\t', header=0)
        df['haplogroup_id'] = df['haplogroup'].apply(lambda x: x[0])

        records = []
        for row in df.itertuples():
            if row.sample == 'NA24385':
                continue
            if row.project == 'REF':
                continue
            sample_info = {
                'sample': row.sample,
                f'is_haplogroup_{row.haplogroup_id}': 1,
                f'is_project_{row.project}': 1
            }
            records.append(sample_info)
            
        features = pd.DataFrame.from_records(records)
        features.fillna(0, inplace=True)
        for col in features.columns:
            if col == 'sample':
                continue
            features[col] = features[col].astype(int)

        features.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule create_complete_feature_table:
    input:
        projects = 'output/quality_factors/project_haplogroup_feature.tsv',
        quast = 'output/eval/assm_stats/SAMPLES.{hifi_type}.{ont_type}.na.chrY.quast-report.tsv',
        qvest = 'output/eval/assembly_qv/SAMPLES.{hifi_type}.{ont_type}.na.wg.yak-qv.tsv',
        reads = 'output/stats/reads/SAMPLES.READS.read-stats.tsv'
    output:
        table = 'output/quality_factors/feature_table.{hifi_type}.{ont_type}.tsv',
    run:
        import pandas as pd

        feat_projects = pd.read_csv(input.projects, sep='\t', header=0)

        feat_quasts = pd.read_csv(
            input.quast, sep='\t', header=0,
            usecols=['sample', 'assembly_length_bp', 'contig_NG50', 'contigs_num', 'largest_contig_bp']
        )
        feat_quasts = feat_quasts.loc[~feat_quasts['sample'].isin(['NA24385']), :].copy()

        feat_table = pd.merge(feat_projects, feat_quasts, on='sample', how='outer')

        feat_qvest = pd.read_csv(
            input.qvest, sep='\t', header=0,
            usecols=['sample', 'qv']
        )
        feat_qvest = feat_qvest.loc[~feat_qvest['sample'].isin(['NA24385']), :].copy()

        feat_table = pd.merge(feat_table, feat_qvest, on='sample', how='outer')

        feat_reads = pd.read_csv(input.reads, sep='\t', header=0)
        feat_reads.sort_values(['sample', 'read_type'], ascending=True, inplace=True)

        sample_read_stats = []
        for row in feat_reads.itertuples():
            if row.sample == 'NA24385':
                continue
            read_type = row.read_type
            assert read_type in [wildcards.hifi_type, wildcards.ont_type]
            if read_type == 'HIFIRW':
                this_sample = dict()
                this_sample['sample'] = row.sample
            this_sample[f'{read_type}_readlen_N50'] = row.read_length_N50_bp
            this_sample[f'{read_type}_readlen_mean'] = row.read_length_mean_bp
            this_sample[f'{read_type}_coverage'] = row.cov_geq_0bp_T2TXYM_linear
            if read_type == 'ONTUL':
                this_sample[f'{read_type}_ul_cov'] = row.cov_geq_100kbp_T2TXYM_linear
                sample_read_stats.append(this_sample)

        sample_read_stats = pd.DataFrame.from_records(sample_read_stats)
        feat_table = pd.merge(feat_table, sample_read_stats, on='sample', how='outer')
        feat_table.sort_values('sample', ascending=True, inplace=True)
        feat_table.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule train_regression_model:
    input:
        table = 'output/quality_factors/feature_table.{hifi_type}.{ont_type}.tsv',
    output:
        model = 'output/quality_factors/model_stats.{hifi_type}.{ont_type}.{target}.tsv',
        data = 'output/quality_factors/model_data.{hifi_type}.{ont_type}.{target}.tsv',
    wildcard_constraints:
        target = '(ctgng50|assmlen|ctgnum)'
    conda:
        '../envs/pyscript.yaml'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = find_script_path('model_assm_quality.py'),
    shell:
        '{params.script_exec} --data-table {input.table} --regression-target {wildcards.target} '
            '--num-cpu {threads} --out-data {output.data} --out-model {output.model}'


rule merge_model_summary_stats:
    input:
        tables = expand(
            'output/quality_factors/model_stats.{{hifi_type}}.{{ont_type}}.{target}.tsv',
            target=['ctgng50', 'assmlen', 'ctgnum']
        )
    output:
        table = 'output/quality_factors/model_stats.{hifi_type}.{ont_type}.all-targets.tsv',
    run:
        import pandas as pd

        model_stats = []
        for table in input.tables:
            df = pd.read_csv(table, header=0, sep='\t')
            model_stats.append(df)

        model_stats = pd.concat(model_stats, axis=0, ignore_index=False)

        model_stats.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK