
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


localrules: collect_all_kmer_stats, collect_all_kmer_differences

rule collect_all_kmer_stats:
    input:
        sample_stats = expand(
            'output/kmer_dump/{sample}.{{hifi_type}}.{{ont_type}}.{{mapq}}.{{chrom}}.k{{kmer}}.meryl-stats.txt',
            sample=SAMPLE_NAMES
        ),
        ref_stats = expand(
            'output/kmer_dump/{reference}.{{chrom}}.k{{kmer}}.meryl-stats.txt',
            reference=['T2T', 'GRCh38'],
        ),
        
    output:
        table = 'output/stats/kmers/SAMPLES.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.counts.tsv'
    run:
        import pandas as pd
        import itertools as itt
        import pathlib as pl

        keywords = ['unique', 'distinct']
        records = []
        for stats_file in itt.chain(input.sample_stats, input.ref_stats):
            assert pl.Path(stats_file).is_file()
            sample = pl.Path(stats_file).stem.split('.')[0]
            sample_record = {
                    'sample': sample
                }
            with open(stats_file, 'r') as table:
                for line in table:
                    if not any(k in line for k in keywords):
                        continue
                    parts = line.strip().split()
                    kmer_type = parts[0]
                    assert kmer_type in keywords
                    type_abundance = int(parts[1])
                    sample_record[kmer_type] = type_abundance
                    if kmer_type == 'distinct':
                        break
            records.append(sample_record)
        df = pd.DataFrame.from_records(records)
        df.sort_values('sample', ascending=True, inplace=True)
        df.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule collect_all_kmer_differences:
    input:
        t2t_hg38 = 'output/eval/kmer_stats/T2T_not-in_GRCh38.{chrom}.k{kmer}.count.txt',
        hg38_t2t = 'output/eval/kmer_stats/GRCh38_not-in_T2T.{chrom}.k{kmer}.count.txt',
        samples = expand(
            'output/eval/kmer_stats/{sample1}_not-in_{sample2}.{{hifi_type}}.{{ont_type}}.{{mapq}}.{{chrom}}.k{{kmer}}.count.txt',
            match_kmer_samples,
            sample1=SAMPLE_NAMES + ['T2T', 'GRCh38'],
            sample2=SAMPLE_NAMES + ['T2T', 'GRCh38'],
        ),
        counts = 'output/stats/kmers/SAMPLES.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.counts.tsv'
    output:
        table = 'output/stats/kmers/SAMPLES.{hifi_type}.{ont_type}.{mapq}.{chrom}.k{kmer}.diffs.tsv'
    run:
        import itertools as itt
        import pathlib as pl
        import pandas as pd

        records = []
        for counts_file in itt.chain([input.t2t_hg38], [input.hg38_t2t], input.samples):
            #assert pl.Path(counts_file).is_file(), f'Not a valid file: {counts_file}'
            trg_sample, qry_sample = pl.Path(counts_file).stem.split('.')[0].split('_not-in_')
            diff_count = int(open(counts_file).read().strip())
            records.append(
                {'target_sample': trg_sample, 'query_sample': qry_sample, 'kmers_target_only': diff_count}
            )
        diffs = pd.DataFrame.from_records(records)
        counts = pd.read_csv(input.counts, sep='\t', header=0,
            dtype={'sample': str, 'unique': int, 'distinct': int}
        )
        diffs = diffs.merge(counts, how='outer', left_on='target_sample', right_on='sample')
        assert pd.notnull(diffs).all(axis=0).all()
        diffs['kmers_target_only_pct'] = (diffs['kmers_target_only'] / diffs['distinct'] * 100).round(3)
        diffs.sort_values(['target_sample', 'query_sample'], ascending=True, inplace=True)
        diffs.drop('sample', axis=1, inplace=True)
        diffs.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


localrules: determine_contiguous_assembly_in_par
rule determine_contiguous_assembly_in_par:
    """
    Because the PAR regions cannot be spanned by a contiguously
    assembled contig, this checks of (i) only a single contig
    has been assembled (ii) to at least 95% of the size in the
    T2T-Y reference
    """
    input:
        t2t = 'references_derived/T2T.chrY-seq-classes.tsv',
        seqclasses = expand(
            'references_derived/{sample}.{{hifi_type}}.{{ont_type}}.na.{chrom}.seqclasses.bed',
            sample=sorted(set(COMPLETE_SAMPLES) - set(['NA24385'])),
            chrom=['chrY']
        )
    output:
        table = 'output/stats/contigs/contig-par.{hifi_type}.{ont_type}.na.chrY.tsv',
    run:
        import pandas as pd
        import pathlib as pl

        t2t = pd.read_csv(input.t2t, sep='\t', header=0)
        t2t['length'] = t2t['end'] - t2t['start']

        sample_stats = []
        for sc_file in input.seqclasses:
            sample = pl.Path(sc_file).stem.split('.')[0]
            df = pd.read_csv(sc_file, sep='\t', header=None, names=['contig', 'start', 'end', 'class_name'])
            df['length'] = df['end'] - df['start']
            stats = {
                'sample': sample,
                'PAR1_contigs': int(df['class_name'].str.contains('PAR1').sum()),
                'PAR2_contigs': int(df['class_name'].str.contains('PAR2').sum())
            }
            for par in ['PAR1', 'PAR2']:
                ref_len = t2t.loc[t2t['name'] == par, 'length'].values[0]
                assm_len = df.loc[df['class_name'] == par, 'length'].values[0]
                stats[f'{par}_assembled_bp'] = assm_len
                pct_assm = round(assm_len / ref_len * 100, 1)
                stats[f'{par}_assembled_pct'] = pct_assm
                stats[f'{par}_is_contiguous'] = 0
                if stats[f'{par}_contigs'] == 1 and pct_assm > 95:
                    stats[f'{par}_is_contiguous'] = 1
            sample_stats.append(stats)

        sample_stats = pd.DataFrame.from_records(sample_stats)
        sample_stats.sort_values(['sample'], inplace=True)
        sample_stats.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


localrules: aggregate_contig_sequence_class_coverage
rule aggregate_contig_sequence_class_coverage:
    input:
        seqclasses = 'references_derived/T2T.chrY-seq-classes.tsv',
        fasta_idx = expand(
            'output/subset_wg/20_extract_contigs/{sample}.{{hifi_type}}.{{ont_type}}.na.chrY.fasta.fai',
            sample=COMPLETE_SAMPLES
        ),
        par_info = 'output/stats/contigs/contig-par.{hifi_type}.{ont_type}.na.chrY.tsv',
    output:
        cov = 'output/stats/contigs/contig-cov.{hifi_type}.{ont_type}.na.chrY.tsv',
        ctg = 'output/stats/contigs/contig-ctg.{hifi_type}.{ont_type}.na.chrY.tsv',
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import numpy as np
        import pathlib as pl

        t2t = pd.read_csv(input.seqclasses, sep='\t', header=0)
        name_idx_map = dict(
            (n,i) for n,i in zip(t2t['name'].values, t2t.index.values)
        )
        num_regions = len(name_idx_map)
        region_names = t2t['name'].values

        ignore_samples = ['HG02666', 'HG01457', 'NA19384', 'NA18989', 'NA24385']

        sample_contigs = []
        contig_covs = []
        contig_ctgs = []

        for idx_file in input.fasta_idx:
            sample = pl.Path(idx_file).stem.split('.')[0]
            if sample in ignore_samples:
                continue
            with open(idx_file, 'r') as faidx:
                for line in faidx:
                    cov = np.zeros(num_regions, dtype=np.int8)
                    ctg = np.zeros(num_regions, dtype=np.int8)

                    contig = line.split()[0]
                    from_class, to_class = contig.split('.')[3].split('-')
                    from_idx = name_idx_map[from_class]
                    to_idx = name_idx_map[to_class]
                    cov[from_idx:to_idx+1] = 1
                    sample_contigs.append(
                        (sample, contig)
                    )
                    contig_covs.append(cov)
                    is_ctg_assm = (to_idx - from_idx) > 1
                    if is_ctg_assm:
                        ctg[from_idx+1:to_idx] = 1
                    contig_ctgs.append(ctg)

        multi_idx = pd.MultiIndex.from_tuples(sample_contigs, names=['sample', 'contig'])
        coverages = pd.DataFrame.from_records(
            contig_covs,
            columns=region_names,
            index=multi_idx
        )
        coverages.sort_index(axis=0, ascending=True, inplace=True)
        coverages.to_csv(output.cov, header=True, index=True, sep='\t')

        contiguity = pd.DataFrame.from_records(
            contig_ctgs,
            columns=region_names,
            index=multi_idx
        )
        # special adaptation: contiguous assembly is defined differently
        # for PAR region, load that info from rule determine_contiguous_assembly_in_par
        par_info = pd.read_csv(input.par_info, sep='\t', header=0)
        for par in ['PAR1', 'PAR2']:
            par_ctg_samples = set(par_info.loc[par_info[f'{par}_is_contiguous'], 'sample'].values)
            select_samples = np.array(
                [s in par_ctg_samples for s in contiguity.index.get_level_values('sample')], dtype=np.bool
            )
            select_contigs = np.array(
                [par in value for value in contiguity.index.get_level_values('contig')], dtype=np.bool
            )
            select_rows = select_samples & select_contigs
            contiguity.loc[select_rows, par] = 1


        # drop sample/contigs w/o contiguous assembly
        contiguity.drop(
            contiguity.index[(contiguity == 0).all(axis=1)],
            axis=0,
            inplace=True,
        )

        # drop PAR1/PAR2 --- no longer needed, see above
        #contiguity.drop(['PAR1', 'PAR2'], axis=1, inplace=True)
        contiguity.sort_index(axis=0, ascending=True, inplace=True)
        contiguity.to_csv(output.ctg, header=True, index=True, sep='\t')
    # END OF RUN BLOCK
