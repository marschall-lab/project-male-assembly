
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
        ref_size = 6017864545, # this is the size for a diploid male assembly (T2T-based) incl. 2 MT
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


rule dump_meryl_chrom_kmer_db:
    input:
        assm = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.fasta',
    output:
        db = directory('output/kmer_dump/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl'),
        stats = 'output/kmer_dump/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl-stats.txt'
    log:
        'log/output/kmer_dump/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl.log'
    benchmark:
        'rsrc/output/kmer_dump/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl.rsrc'
    wildcard_constraints:
        sample = SAMPLE_NAME_CONSTRAINT
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, input, attempt: 2048 * attempt,
        mem_gb = lambda wildcards, input, attempt: 2 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    params:
        kmer_size = lambda wildcards: int(wildcards.kmer)
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_gb} output {output.db} {input.assm} &> {log}'
            ' && '
        'meryl statistics {output.db} > {output.stats}'


rule dump_meryl_ref_chrom_kmer_db:
    input:
        assm = 'references_derived/{reference}_chrY.fasta',
    output:
        db = directory('output/kmer_dump/{reference}.chrY.k{kmer}.meryl'),
        stats = 'output/kmer_dump/{reference}.chrY.k{kmer}.meryl-stats.txt'
    log:
        'log/output/kmer_dump/{reference}.chrY.k{kmer}.meryl.log'
    benchmark:
        'rsrc/output/kmer_dump/{reference}.chrY.k{kmer}.meryl.rsrc'
    wildcard_constraints:
        reference = '(T2T|GRCh38)'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, input, attempt: 2048 * attempt,
        mem_gb = lambda wildcards, input, attempt: 2 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    params:
        kmer_size = lambda wildcards: int(wildcards.kmer)
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_gb} output {output.db} {input.assm} &> {log}'
            ' && '
        'meryl statistics {output.db} > {output.stats}'


rule difference_meryl_sample_kmer_dbs:
    input:
        s1 = 'output/kmer_dump/{sample1}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl',
        s2 = 'output/kmer_dump/{sample2}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl'
    output:
        s1_only = 'output/eval/kmer_stats/{sample1}_not-in_{sample2}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.count.txt',
        s2_only = 'output/eval/kmer_stats/{sample2}_not-in_{sample1}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.count.txt',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample1 = SAMPLE_NAME_CONSTRAINT,
        sample2 = SAMPLE_NAME_CONSTRAINT
    resources:
        mem_mb = lambda wildcards, input, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    shell:
        'meryl print [difference {input.s1} {input.s2}] | wc -l > {output.s1_only}'
            ' && '
        'meryl print [difference {input.s2} {input.s1}] | wc -l > {output.s2_only}'


rule difference_meryl_ref_sample_kmer_dbs:
    input:
        smp = 'output/kmer_dump/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.meryl',
        ref = 'output/kmer_dump/{reference}.{chrom}.k{kmer}.meryl'
    output:
        smp_only = 'output/eval/kmer_stats/{sample}_not-in_{reference}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.count.txt',
        ref_only = 'output/eval/kmer_stats/{reference}_not-in_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.count.txt',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        reference = '(GRCh38|T2T)',
        sample = SAMPLE_NAME_CONSTRAINT
    resources:
        mem_mb = lambda wildcards, input, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    shell:
        'meryl print [difference {input.smp} {input.ref}] | wc -l > {output.smp_only}'
            ' && '
        'meryl print [difference {input.ref} {input.smp}] | wc -l > {output.ref_only}'


rule difference_meryl_ref_ref_kmer_dbs:
    input:
        ref1 = 'output/kmer_dump/{ref1}.{chrom}.k{kmer}.meryl',
        ref2 = 'output/kmer_dump/{ref2}.{chrom}.k{kmer}.meryl'
    output:
        ref1_only = 'output/eval/kmer_stats/{ref1}_not-in_{ref2}.{chrom}.k{kmer}.count.txt',
        ref2_only = 'output/eval/kmer_stats/{ref2}_not-in_{ref1}.{chrom}.k{kmer}.count.txt',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        ref1 = '(GRCh38|T2T)',
        ref2 = '(GRCh38|T2T)',
        chrom = 'chrY'
    resources:
        mem_mb = lambda wildcards, input, attempt: 2048 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00'
    shell:
        'meryl print [difference {input.ref1} {input.ref2}] | wc -l > {output.ref1_only}'
            ' && '
        'meryl print [difference {input.ref2} {input.ref1}] | wc -l > {output.ref2_only}'


def match_kmer_samples(*args, **kwargs):

    samples1 = [s for _, s in args[0]]
    samples2 = [s for _, s in args[1]]

    refs = ['T2T', 'GRCh38']

    pairs = []
    for s1 in samples1:
        for s2 in samples2:
            if s1 == s2:
                continue
            if s1 in refs and s2 in refs:
                continue
            pairs.append(
                {'sample1': s1, 'sample2': s2}
            )
    return pairs

localrules: collect_all_kmer_stats
rule collect_all_kmer_stats:
    input:
        t2t_hg38 = 'output/eval/kmer_stats/T2T_not-in_GRCh38.chrY.k31.count.txt',
        hg38_t2t = 'output/eval/kmer_stats/GRCh38_not-in_T2T.chrY.k31.count.txt',
        samples = expand(
            'output/eval/kmer_stats/{sample1}_not-in_{sample2}.HIFIRW.ONTUL.na.chrY.k31.count.txt',
            match_kmer_samples,
            sample1=SAMPLE_NAMES + ['T2T', 'GRCh38'],
            sample2=SAMPLE_NAMES + ['T2T', 'GRCh38'],
        )
    output:
        touch('output/eval/kmer_stats/all-vs-all.ok')
