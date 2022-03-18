
"""
Only motifs that are used for Y contig identification are thresholded;
for all other motifs, all chrY hits are retained

- DYZ1 - score >2500
- DYZ2 - score >1700
- DYZ18 - score >2100
- DYZ3-sec_Ycentro - score >1700
"""

RUNTIME_MOTIF_FACTOR = {
    'DYZ2_Yq': 10,
}

CPU_MOTIF_FACTOR = {
    'DYZ2_Yq': config['num_cpu_high'],
}

BONUS_MOTIF_FACTOR = {
    'DYZ2_Yq': 100,
}

rule hmmer_motif_search:
    input:
        assm = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        qry = 'references_derived/{motif}.fasta'
    output:
        txt = 'output/motif_search/00_detection/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.txt',
        table = 'output/motif_search/00_detection/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
    log:
        'log/output/motif_search/00_detection/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.log',
    benchmark:
        'rsrc/output/motif_search/00_detection/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.rsrc',
#    singularity:
#        'hmmer.sif'
    conda:
        '../envs/biotools.yaml'
    threads: lambda wildcards: CPU_MOTIF_FACTOR.get(wildcards.motif, config['num_cpu_medium'])
    resources:
        mem_mb = lambda wildcards, attempt: 24576 + 24576 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt*RUNTIME_MOTIF_FACTOR.get(wildcards.motif, 1):02}:59:00',
        bonus = lambda wildcards, attempt: BONUS_MOTIF_FACTOR.get(wildcards.motif, 0)
    params:
        evalue = '1.60E-150'
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

SCORE_THRESHOLDS_MOTIF = {
    'DYZ1_Yq': 2500,
    'DYZ18_Yq': 2100,
    'DYZ2_Yq': 1700,
    'DYZ3-sec_Ycentro': 1700
}

rule normalize_motif_hits:
    input:
        table = 'output/motif_search/00_detection/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
        motif = 'references_derived/{motif}.fasta',
    output:
        table = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
        bed_all = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.bed',
        bed_hiq = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm-hiq.bed',
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        min_score_t = lambda wildcards: SCORE_THRESHOLDS_MOTIF.get(wildcards.motif, 0)
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
    input:
        table = 'output/motif_search/10_norm/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.norm.tsv',
    output:
        table = 'output/motif_search/20_target_agg/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.agg-trg.tsv',
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
