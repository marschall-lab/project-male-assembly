
"""
DYZ1 (Y-specific), DYZ18 (Y-specific, but found only at the AMPL7 and Yqhet boundary) and DYZ2 (found on Y and elsewhere). 

I filter the nhmmer_out_parse.txt based on the 'score' column (column 14):
- DYZ1 - score >2500
- DYZ2 - score >1700
- DYZ18 - score >2100
"""

rule hmmer_motif_search:
    input:
        assm = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
        qry = 'references_derived/{motif}.fasta'
    output:
        txt = 'output/motif_search/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.txt',
        table = 'output/motif_search/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.table.txt',
    log:
        'log/output/motif_search/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.log',
    benchmark:
        'rsrc/output/motif_search/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.{motif}.hmmer.rsrc',
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 49152 + 24576 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt**3:02}:59:00',
    params:
        evalue = '1.60E-150'
    shell:
        'nhmmer --cpu {threads} -o {output.txt} --tblout {output.table} -E {params.evalue} {input.qry} {input.assm} &> {log}'
