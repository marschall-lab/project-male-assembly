

rule convert_contig_paf_to_bed:
    input:
        paf = 'output/alignments/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    output:
        bed = 'output/eval/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed'
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    shell:
        'zgrep -F "tp:A:P" {input.paf} | cut -f 1,5,6,8,9,12 | '
            'awk \'BEGIN{{OFS = "\\t"}} {{print $3,$4,$5,$1,$6,$2}}\' | '
            'sort -V -k1 -k2n,3n > {output}'


rule intersect_contig_align_chry_seq_classes:
    input:
        align_bed = 'output/eval/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed',
        ref_bed = 'references_derived/{reference}.bed'
    output:
        tsv = 'output/eval/contigs-to-ref/{sub_folder}/{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.ref-cov.tsv'
    wildcard_constraints:
        sub_folder = '(00_raw|10_renamed)'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'bedtools intersect -wao -a {input.ref_bed} -b {input.align_bed} > {output.tsv}'


rule determine_majority_haplotype_hifiasm:
    input:
        paf = 'output/alignments/contigs-to-contigs/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_hifiasm.paf.gz',
    output:
        table = 'output/eval/contigs-to-contigs/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_hifiasm.stats.tsv',
    run:
        import pandas as pd
        import numpy as np

        paf_columns = ['qname', 'qlen', 'qstart', 'qend', 'strand',
            'tname', 'tlen', 'tstart', 'tend', 'nmatch', 'blen', 'mapq',
            'tag_nm', 'tag_ms', 'tag_as', 'tag_nn', 'tag_tp' ]

        df = pd.read_csv(input.paf, sep='\t', header=None, names=paf_columns, usecols=range(len(paf_columns)))
        df['tag_tp'] = df['tag_tp'].str.lower()

        df['hap'] = df['tname'].str.slice(0,2)  # specific to tig-naming in hifiasm
        # fix for HG01890 - one unaligned contig
        rename_haps = {
            'h1': 'h1',
            'h2': 'h2'
        }
        df['hap'] = df['hap'].apply(lambda x: rename_haps.get(x, 'un'))
        assert (df['hap'].isin(['h1', 'h2', 'un'])).all()

        records = []  # collect records for flat table
        for hap, aligns in df.groupby('hap'):
            records.extend([
                (hap, 'all', 'alignments_total', aligns.shape[0]),
                (hap, 'all', 'queries_total', aligns['qname'].nunique())
            ])
            if hap == 'un':
                continue
            
            is_primary = aligns['tag_tp'].isin(['tp:a:p', 'tp:a:i'])
            aln = aligns.loc[is_primary, :].copy()
            records.extend([
                (hap, 'all', 'alignments_not_2nd', aln.shape[0]),
                (hap, 'all', 'queries_not_2nd', aln['qname'].nunique())
            ])
            
            aln['qalign_len'] = aln['qend'] - aln['qstart']    
            avg_mapq = np.average(aln['mapq'], weights=aln['qalign_len'])
            records.append(
                (hap, 'all', 'wtavg_mapq', avg_mapq)
            )

            is_hq = aln['mapq'] == 60
            aln = aln.loc[is_hq, :].copy()
            records.extend([
                (hap, 'all', 'alignments_is_hq', aln.shape[0]),
                (hap, 'all', 'queries_is_hq', aln['qname'].nunique())
            ])

            aligned_contigs = set()
            for contig, ctg_aln in aln.groupby('qname'):
                aligned_contigs.add(contig)

                contig_size = ctg_aln.at[ctg_aln.index[0], 'qlen']
                sum_aligned = ctg_aln['qalign_len'].sum()
                hasm_contigs = ctg_aln['tname'].nunique()
                hasm_contig_lens = ctg_aln[['tname', 'tlen']].drop_duplicates(inplace=False)['tlen'].sum()
                records.extend([
                    (hap, contig, 'query_size', contig_size),
                    (hap, contig, 'query_aligned_bp', sum_aligned),
                    (hap, contig, 'query_aligned_pct', round(sum_aligned / contig_size * 100, 2)),
                    (hap, contig, 'target_contigs', hasm_contigs),
                    (hap, contig, 'target_length_total', hasm_contig_lens),
                    (hap, contig, 'alignment_spread', hasm_contig_lens / contig_size),            
                ])

            for contig, ctg_aln in aligns.loc[~aligns['qname'].isin(aligned_contigs), :].groupby('qname'):
                contig_size = ctg_aln.at[ctg_aln.index[0], 'qlen']
                records.extend([
                    (hap, contig, 'query_size', contig_size),
                    (hap, contig, 'query_aligned_bp', 0),
                    (hap, contig, 'query_aligned_pct', 0),
                    (hap, contig, 'target_contigs', 0),
                    (hap, contig, 'target_length_total', 0),
                    (hap, contig, 'alignment_spread', 0),            
                ])

        max_hap = None
        max_avg = 0
        for hap in ['h1', 'h2']:
            hap_aln_pct = [t[3] for t in records if t[0] == 'h1' and t[1] != 'all' and t[2] == 'query_aligned_pct']
            hap_aln_wt = [t[3] for t in records if t[0] == 'h1' and t[1] != 'all' and t[2] == 'query_size']
            wtavg = np.average(hap_aln_pct, weights=hap_aln_wt)
            records.append((hap, 'all', 'wtavg_qalign_pct', wtavg))
            if wtavg > max_avg:
                max_hap = hap
                max_avg = wtavg
        records.append((max_hap, 'all', 'select_majority_hap', int(max_hap.strip('h'))))

        summary = pd.DataFrame.from_records(
            records,
            columns=['hap', 'contig', 'statistic', 'value']
        )
        summary.sort_values(['contig', 'hap', 'statistic'], ascending=True, inplace=True)
        summary.to_csv(output.table, sep='\t', header=True, index=False)
    # END OF RUN BLOCK