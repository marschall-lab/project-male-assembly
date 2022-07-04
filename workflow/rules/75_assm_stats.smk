
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

        keep_records = dict(
            [(1, 'sample'), (2, 'contigs_num'), (7, 'contigs_geq50kbp_num'),
            (8, 'assembly_len_bp'), (15, 'largest_contig'), (17, 'reference_len_bp'), (18, 'GC_pct'),
            (19, 'contig_N50'), (20, 'contig_NG50'), (21, 'contig_N90'), (22, 'contig_auN'),
            (23, 'contig_auNG'), (24, 'contig_L50'), (25, 'contig_LG50'), (26, 'contig_L90')]
        )

        rows = []
        for file_path in input.reports:
            with open(file_path, 'r') as report:
                row = dict()
                for ln, line in enumerate(report, start=1):
                    try:
                        column_name = keep_records[ln]
                        column_value = line.strip().split()[-1]
                        if column_name == 'sample':
                            pass
                        elif column_name in ['GC_pct', 'contig_auN', 'contig_auNG']:
                            column_value = float(column_value)
                        else:
                            column_value = int(column_value)
                    except KeyError:
                        continue
                    row[column_name] = column_value
            rows.append(row)

        df = pd.DataFrame.from_records(rows)
        df.sort_values('sample', ascending=True, inplace=True)

        df.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK
