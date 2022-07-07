

RUNTIME_MOTIF_FACTOR = {
    'DYZ2_Yq': 10,
    'DYZ3-sec_Ycentro': 2,
    'TSPY': 35
}

CPU_MOTIF_FACTOR = {
    'DYZ2_Yq': config['num_cpu_high'],
    'TSPY': config['num_cpu_high'] + config['num_cpu_medium']
}

MEMORY_MOTIF_FACTOR = {
    'TSPY': 15
}


rule hmmer_motif_search:
    input:
        assm = select_whole_genome_assembly,
        qry = 'references_derived/{motif}.fasta'
    output:
        txt = 'output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.txt',
        table = 'output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
    log:
        'log/output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.log',
    benchmark:
        'rsrc/output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.rsrc',
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    singularity:
        'hmmer.sif'
#    conda:
#        '../envs/biotools.yaml'
    threads: lambda wildcards: CPU_MOTIF_FACTOR.get(wildcards.motif, config['num_cpu_medium'])
    resources:
        mem_mb = lambda wildcards, attempt: (24576 + 24576 * attempt) * MEMORY_MOTIF_FACTOR.get(wildcards.motif, 1),
        walltime = lambda wildcards, attempt: f'{attempt*RUNTIME_MOTIF_FACTOR.get(wildcards.motif, 1):02}:59:00',
    params:
        evalue = lambda wildcards: config['hmmer_evalue_cutoff'][wildcards.motif]
    shell:
        'nhmmer --cpu {threads} --dna -o {output.txt} --tblout {output.table} -E {params.evalue} {input.qry} {input.assm} &> {log}'


rule hmmer_reference_motif_search:
    input:
        assm = 'references_derived/{sample}_chrY.fasta',
        qry = 'references_derived/{motif}.fasta'
    output:
        txt = 'output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.txt',
        table = 'output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
    log:
        'log/output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.log',
    benchmark:
        'rsrc/output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.rsrc',
    wildcard_constraints:
        sub_folder = '20_refseq'
    singularity:
        'hmmer.sif'
#    conda:
#        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
    params:
        evalue = lambda wildcards: config['hmmer_evalue_cutoff'][wildcards.motif]
    shell:
        'nhmmer --cpu {threads} --dna -o {output.txt} --tblout {output.table} -E {params.evalue} {input.qry} {input.assm} &> {log}'


HMMER_TABLE_COLUMNS = [
    ('target', True),
    ('target_accession', False),
    ('query', True),
    ('query_accession', False),
    ('query_hit_start', True),
    ('query_hit_end', True),
    ('target_hit_start', True),
    ('target_hit_end', True),
    ('target_env_start', True),
    ('target_env_end', True),
    ('target_length', True),
    ('target_strand', True),
    ('evalue', True,),
    ('bit_score', True),
    ('bias', True),
    ('description', False)
]

HMMER_TABLE_NAMES = [t[0] for t in HMMER_TABLE_COLUMNS]
HMMER_TABLE_USE_COLS = [t[0] for t in HMMER_TABLE_COLUMNS if t[1]]

rule normalize_motif_hits:
    input:
        table = 'output/motif_search/00_detection/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
        motif = 'references_derived/{motif}.fasta',
    output:
        table = 'output/motif_search/10_norm/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
        bed_all = 'output/motif_search/10_norm/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.bed',
        bed_hiq = 'output/motif_search/10_norm/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm-hiq.bed',
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed|20_refseq)'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        min_score_t = lambda wildcards: config['hmmer_score_threshold'].get(wildcards.motif, 0)
    run:
        import pandas as pd

        query_length = get_fasta_seq_length(input.motif)[wildcards.motif]

        df = pd.read_csv(
            input.table,
            delimiter=r"\s+",
            header=None,
            skip_blank_lines=True,
            comment='#',
            names=HMMER_TABLE_NAMES,
            usecols=HMMER_TABLE_USE_COLS
        )
        df['query_length'] = query_length
        df['query_hit_pct'] = (df['query_hit_end'] - (df['query_hit_start'] - 1)) / df['query_length'] * 100
        df['query_hit_pct'] = df['query_hit_pct'].round(2)
        df['hit_hiq'] = 0
        df.loc[df['bit_score'] > params.min_score_t, 'hit_hiq'] = 1

        # for hits on the revcomp strand, switch start and end to create valid BED output
        revcomp_select = df['target_strand'] == '-'
        forward_cols = ['target_hit_start', 'target_hit_end', 'target_env_start', 'target_env_end']
        reverse_cols = ['target_hit_end', 'target_hit_start', 'target_env_end', 'target_env_start']
        df.loc[revcomp_select, forward_cols] = df.loc[revcomp_select, reverse_cols].values
        assert (df['target_hit_start'] < df['target_hit_end']).all()

        df.to_csv(output.table, sep='\t', header=True, index=False)
        bed_columns = ['target', 'target_hit_start', 'target_hit_end', 'query', 'bit_score', 'target_strand', 'hit_hiq']
        df = df[bed_columns]
        df['target_hit_start'] -= 1
        df.sort_values(['target', 'target_hit_start', 'target_hit_end'], inplace=True)
        df.to_csv(output.bed_all, sep='\t', header=False, index=False, columns=bed_columns[:-1])

        df = df.loc[df['hit_hiq'] > 0, :]
        df.to_csv(output.bed_hiq, sep='\t', header=False, index=False, columns=bed_columns[:-1])
    # END OF RUN BLOCK


rule aggregate_motif_hits_by_target:
    """
    Change 2022-05-10
    Rare case: sample HG02572 does not have any high-quality matches for DYZ18, hence
    the original code produced an error. Change to always have a consistent output table.
    """
    input:
        table = 'output/motif_search/10_norm/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
    output:
        table = 'output/motif_search/20_target_agg/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.agg-trg.tsv',
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    run:
        import pandas as pd
        import collections as col

        df = pd.read_csv(input.table, sep='\t', header=0)

        records = []
        for trg, hits in df.groupby('target'):
            record = col.OrderedDict({
                'target': trg,
                'query': wildcards.motif,
                'target_length': hits.at[hits.index[0], 'target_length'],
                'num_hits_total': hits.shape[0],
                'num_hits_hiq': hits.loc[hits['hit_hiq'] > 0, :].shape[0],
                'num_hits_loq': hits.loc[hits['hit_hiq'] < 1, :].shape[0],
            })
            if record['num_hits_hiq'] == 0:
                # 2022-05-10 add empty fields for output table
                record['pct_hits_hiq'] = 0
                record['num_bp_hiq'] = 0
                record['pct_bp_hiq'] = 0
                record['mean_pct_len_hiq'] = 0
                record['median_pct_len_hiq'] = 0
                record['mean_score_hiq'] = 0
                record['median_score_hiq'] = 0
                records.append(record)
                continue
            record['pct_hits_hiq'] = round(record['num_hits_hiq'] / record['num_hits_total'] * 100, 2)
            hiq_subset = hits.loc[hits['hit_hiq'] > 0, :].copy()
            record['num_bp_hiq'] = (hiq_subset['target_hit_end'] - (hiq_subset['target_hit_start'] - 1)).sum()
            record['pct_bp_hiq'] = round(record['num_bp_hiq'] / record['target_length'] * 100, 2)
            record['mean_pct_len_hiq'] = round(hiq_subset['query_hit_pct'].mean(), 2)
            record['median_pct_len_hiq'] = round(hiq_subset['query_hit_pct'].median(), 2)
            record['mean_score_hiq'] = round(hiq_subset['bit_score'].mean(), 0)
            record['median_score_hiq'] = round(hiq_subset['bit_score'].median(), 0)
            records.append(record)

        agg = pd.DataFrame.from_records(records).fillna(0, inplace=False)
        agg['num_bp_hiq'] = agg['num_bp_hiq'].astype(int)
        agg.sort_values('target', inplace=True)
        agg.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule repmsk_chry_contigs:
    """
    With renamed chrY contigs, RepeatMasker fails b/c of...
    $ Fasta file contains a sequence identifier which is too long ( max id length = 50 )
    so create a temp copy with original names, and rename output later
    """
    input:
        fasta = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta',
        rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.nto-map.sed'
    output:
        tmp_fasta = temp('output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta'),
        repmask_output = multiext(
            'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY',
            '.fasta.cat.gz', '.fasta.masked', '.fasta.out', '.fasta.tbl'
        )
    log:
        'log/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.log',
    benchmark:
        'rsrc/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.rsrc',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample = SAMPLE_NAME_CONSTRAINT
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
    params:
        out_dir = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY'
    shell:
        'sed -f {input.rename} {input.fasta} > {output.tmp_fasta} && '
        'RepeatMasker -pa {threads} -s -dir {params.out_dir} -species human {output.tmp_fasta} &> {log}'


rule repmsk_ref_chry_contigs:
    """
    Same as for HMMER above.
    For reference chrY, no renaming needed, that part is just kept here
    to not break anything else
    """
    input:
        fasta = 'references_derived/{sample}_chrY.fasta',
        rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.nto-map.sed'
    output:
        tmp_fasta = temp('output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta'),
        repmask_output = multiext(
            'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY',
            '.fasta.cat.gz', '.fasta.masked', '.fasta.out', '.fasta.tbl'
        )
    log:
        'log/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.log',
    benchmark:
        'rsrc/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.rsrc',
    conda:
        '../envs/biotools.yaml'
    wildcard_constraints:
        sample = '(T2T|GRCh38)'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
    params:
        out_dir = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY'
    shell:
        'sed -f {input.rename} {input.fasta} > {output.tmp_fasta} && '
        'RepeatMasker -pa {threads} -s -dir {params.out_dir} -species human {output.tmp_fasta} &> {log}'


rule collect_repeatmasker_output:
    input:
        report = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.tbl',
        masked = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.masked',
        aln = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.cat.gz',
        table = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.out',
        rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.sed',
    output:
        fasta_rn = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.rm-mask.fasta',
        result_tar = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.rm-out.tar.gz',
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        tar_dir = lambda wildcards, input: pathlib.Path(input.report).parent,
        tar_report = lambda wildcards, input: pathlib.Path(input.report).name,
        tar_masked = lambda wildcards, input: pathlib.Path(input.masked).name,
        tar_aln = lambda wildcards, input: pathlib.Path(input.aln).name,
        tar_table = lambda wildcards, input: pathlib.Path(input.table).name,
    shell:
        'seqtk seq -C -A {input.masked} | sed -f {input.rename} > {output.fasta_rn}'
            ' && '
        'tar czf {output.result_tar} -C {params.tar_dir} ./{params.tar_report} ./{params.tar_masked} '
            './{params.tar_aln} ./{params.tar_table}'


rule normalize_repeatmasker_table:
    input:
        table = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta.out',
        rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.otn-map.json',
    output:
        tsv = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.na.chrY.matches.tsv',
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import json

        names_explained = [
            (0, 'sw_score', 'Smith-Waterman score of the match'),
            (1, 'pct_subs', 'Pct. substitutions in matching region compared to the consensus'),
            (2, 'pct_del', 'Pct. of bases opposite a gap in the query sequence (deleted bp)'),
            (3, 'pct_ins', 'Pct. of bases opposite a gap in the repeat consensus (inserted bp)'),
            (4, 'query_name', 'Query sequence name'),
            (5, 'query_start', 'Starting position of match in query sequence'),
            (6, 'query_end', 'Ending position of match in query sequence'),
            (7, 'query_after_bp', 'bp in query sequence past the ending position of match'),
            (8, 'is_complement_match', 'True: match is with the Complement of the consensus sequence in the database'),
            (9, 'repeat_name', 'Name of the matching interspersed repeat'),
            (10, 'repeat_class', 'Class of the repeat'),
            (11, 'repeat_start', 'Starting position of match in database sequence (using top-strand numbering)'),
            (12, 'repeat_end', 'Ending position of match in database sequence'),
            (13, 'repeat_before_bp', 'bp in (complement of) the repeat consensus sequence prior to beginning of the match'),
            (14, 'match_ID', 'Consecutive ID of match'),
            (15, 'is_partly_included', 'True: there is a higher-scoring match whose domain partly (<80%) includes the domain of this match')
        ]

        df = pd.read_csv(
            input.table,
            header=None,
            index_col=False,
            delimiter=r"\s+",
            skip_blank_lines=True,
            skiprows=[0,1],
            comment='#',
            names=[n[1] for n in names_explained]
        )

        name_map = json.loads(open(input.rename, 'r').read())

        for c in ['query_after_bp', 'repeat_start', 'repeat_before_bp']:
            df[c] = df[c].apply(lambda x: int(x.strip('()')))
            df[c] = df[c].astype(int)

        df['is_complement_match'] = df['is_complement_match'].apply(lambda x: True if x == 'C' else False)
        df['is_complement_match'] = df['is_complement_match'].astype(bool)
        df['is_partly_included'] = df['is_partly_included'].apply(lambda x: True if x == '*' else False)
        df['is_partly_included'] = df['is_partly_included'].astype(bool)
        df['query_name'].replace(name_map, inplace=True)
        df.sort_values(['query_name', 'query_start', 'sw_score'], ascending=[True, True, False], inplace=True)

        with open(output.tsv, 'w') as table:
            for pos, head, comment in names_explained:
                _ = table.write(f'##.{pos} - {head} - {comment}\n')
            df.to_csv(table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK
