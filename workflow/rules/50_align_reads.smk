
# rule align_subset_to_chrom:
#     """
#     This rule only aligns read subsets
#     to a chromosome assembly (Y, potentially, X and XY),
#     but for all possible combinations. The most important
#     one is: the assembly reads, i.e. HIFIAF (by default)
#     and ONTUL
#     """
#     input:
#         ctg = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly.fasta',
#         reads = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{other_reads}.{chrom}-reads.{mapq}.fasta.gz',
#     output:
        