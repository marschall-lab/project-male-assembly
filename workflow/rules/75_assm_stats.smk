
localrules: aggregate_quast_reports

rule compute_wg_assembly_stats:
    input:
        assm = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta'
    output:
        report = 'output/eval/assm_stats/quast/{sample}.{hifi_type}.{ont_type}.na.wg/report.tsv',
        transposed = 'output/eval/assm_stats/quast/{sample}.{hifi_type}.{ont_type}.na.wg/transposed_report.tsv',
    singularity:
        f'{config["container_store"]}/{config["container"]["quast"]}'
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 6144 * attempt
    params:
        ref_size = 6017864545, # this is the size for a diploid male assembly (T2T-based) incl. 1 MT
        out_dir = lambda wildcards, output: pathlib.Path(output.report).parent
    shell:
        'quast --no-plots --no-html --no-icarus --min-contig 1 '
            '--output-dir {params.out_dir} --labels {wildcards.sample} '
            '--est-ref-size {params.ref_size} {input.assm}'


rule compute_chrom_assembly_stats:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.fasta'
    output:
        report = 'output/eval/assm_stats/quast/{sample}.{hifi_type}.{ont_type}.na.{chrom}/report.tsv',
        transposed = 'output/eval/assm_stats/quast/{sample}.{hifi_type}.{ont_type}.na.{chrom}/transposed_report.tsv',
    singularity:
        f'{config["container_store"]}/{config["container"]["quast"]}'
    wildcard_constraints:
        chrom = '(chrY|chrX)'
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        ref_size = lambda wildcards: {'chrY': 62460029, 'chrX': 154259566}[wildcards.chrom],
        out_dir = lambda wildcards, output: pathlib.Path(output.report).parent
    shell:
        'quast --no-plots --no-html --no-icarus --min-contig 1 '
            '--output-dir {params.out_dir} --labels {wildcards.sample} '
            '--est-ref-size {params.ref_size} {input.assm}'


rule aggregate_quast_reports:
    """
    Turns out, there is a bug in Quast that
    silently omits certain statistics (so far
    observed: NG90, LG90) for some samples for
    unclear reasons. This leads to a change in
    line numbering, requiring much more bloated
    code to extract the relevant stats.
    """
    input:
        reports = expand(
            'output/eval/assm_stats/quast/{sample}.{{hifi_type}}.{{ont_type}}.na.{{chrom}}/report.tsv',
            sample=SAMPLE_NAMES
        )
    output:
        table = 'output/eval/assm_stats/SAMPLES.{hifi_type}.{ont_type}.na.{chrom}.quast-report.tsv'
    wildcard_constraints:
        chrom = '(wg|chrY|chrX)'
    run:
        import pandas as pd
        import collections as col
        
        keep_records = [
            ('Assembly', 'sample', str, None), ('# contigs (>= 0 bp)', 'contigs_num', int, None),
            ('# contigs (>= 50000 bp)', 'contigs_geq50kbp_num', int, None),
            ('Total length (>= 0 bp)', 'assembly_length_bp', int, None), ('Largest contig', 'largest_contig_bp', int, None),
            ('Estimated reference length', 'reference_length_bp', int, None), ('GC (%)', 'GC_pct', float, None),
            ('auN', 'contig_auN', float, -1), ('auNG', 'contig_auNG', float, -1)
        ]
        for contig_stat in ['N50', 'NG50', 'N90', 'NG90', 'L50', 'L90', 'LG50', 'LG90']:
            keep_records.append((contig_stat, f'contig_{contig_stat}', int, -1))
        keep_records.append(None)

        rows = []
        for file_path in input.reports:
            with open(file_path, 'r') as report:
                this_records = col.deque(keep_records)
                row = dict()
                for ln, line in enumerate(report, start=1):
                    while 1:
                        keep_line = this_records.popleft()
                        if keep_line is None:
                            this_records.append(None)
                            break
                        if line.startswith(keep_line[0]):
                            row_name = keep_line[1]
                            value_type = keep_line[2]
                            row_value = value_type(line.strip().split()[-1])
                            row[row_name] = row_value
                            continue
                        else:
                            this_records.append(keep_line)
                            continue
                # add defaults if stats are missing
                for missing_stat in this_records:
                    if missing_stat is None:
                        continue
                    default_value = missing_stat[3]
                    if default_value is None:
                        raise ValueError(f'Missing statistic in QUAST report w/o default: {file_path} / {missing_stat}')
                    row[missing_stat[1]] = default_value

            rows.append(row)

        df = pd.DataFrame.from_records(rows)
        df['assembly_type'] = wildcards.chrom
        df['is_E2E_ref'] = 0
        select_e2e_ref = df['reference_length_bp'] * 0.99 <= df['contig_NG50']
        df.loc[select_e2e_ref, 'is_E2E_ref'] = 1
        df['is_E2E_assm'] = 0
        select_e2e_assm = df['assembly_length_bp'] * 0.99 <= df['contig_N50']
        df.loc[select_e2e_assm, 'is_E2E_assm']  = 1
        df.sort_values('sample', ascending=True, inplace=True)

        df.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK
