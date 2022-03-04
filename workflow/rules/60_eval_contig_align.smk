

rule convert_contig_paf_to_bed:
    input:
        paf = 'output/alignments/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.paf.gz'
    output:
        bed = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed'
    shell:
        'zgrep -F "tp:A:P" {input.paf} | cut -f 1,5,6,8,9,12 | '
            'awk \'BEGIN{{OFS = "\\t"}} {{print $3,$4,$5,$1,$6,$2}}\' | '
            'sort -V -k1 -k2n,3n > {output}'


rule intersect_contig_align_chry_seq_classes:
    input:
        align_bed = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.bed',
        ref_bed = 'references_derived/{reference}.bed'
    output:
        tsv = 'output/eval/contigs-to-ref/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}_aln-to_{reference}.ref-cov.tsv'
    conda:
        '../envs/biotools.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        'bedtools intersect -wao -a {input.ref_bed} -b {input.align_bed} > {output.tsv}'