
########################################################################################
# Update 2022-05-25
# VerityMap (v2.1.1-alpha) can only be set up inside a Singularity container,
# because of an include error that appears when pip install is run
# in an (otherwise identical) Conda environment (see gh#18). However, the
# containerized version fails when processig wg assemblies w/o producing an
# informative error message, which makes debugging impossible. Hence,
# the ugly setup strategy below: build the environment in the container,
# but run setup outside of it to support code inspection (changes) on-the-fly.
# This creates the veritymap script in a .local cache folder, and other abominations...
# obvious TODO: replace that when errors are fixed
########################################################################################

localrules: clone_veritymap_repo, build_veritymap, clean_veritymap_build, merge_error_annotations


rule clone_veritymap_repo:
    """
    2022-07-04
    Checkout target is bugfix branch
    """
    output:
        'repos/veritymap.git.ok'
    envmodules:
        'git/2.24.0'
    params:
        veritymap_repo = config['veritymap_repo']
    shell:
        'mkdir -p repos/ && cd repos/'
            ' && '
        'git clone {params.veritymap_repo}'
            ' && '
        'cd VerityMap'
            ' && '
        'git checkout 8d241f4'
            ' && '
        'cd ../../ && touch {output}'


rule clean_veritymap_build:
    input:
        'repos/veritymap.git.ok'
    output:
        'repos/veritymap.clean.ok'
    shell:
        'rm -rfd .cache/pip'
            ' && '
        'rm -rfd .local/lib && rm -rfd .local/bin'
            ' && '
        'rm -rfd repos/VerityMap/.cache'
            ' && '
        'rm -rfd repos/VerityMap/VerityMap.egg-info'
            ' && '
        'rm -rfd repos/VerityMap/build'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/build'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/__pycache__'
            ' && '
        'rm -rfd repos/VerityMap/veritymap/py_src/__pycache__'
            ' && '
        'touch {output}'


rule build_veritymap:
    input:
        'repos/veritymap.git.ok',
        'repos/veritymap.clean.ok'
    output:
        'repos/veritymap.build.ok'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap_env"]}'
    shell:
        'pip install --no-warn-script-location --no-deps repos/VerityMap/'
            ' && '
        'touch {output}'


rule run_chrom_veritymap:
    """
    NB: important to execute the main.py at the right
    location due to the setup horror explained above.
    Since main is messing with sys.path, VerityMap
    can be properly executed irrespective of the install/setup
    into "." locations by pip.
    """
    input:
        build_ok = 'repos/veritymap.build.ok',
        reads = 'output/subset_wg/45_extract_reads/{sample}.{other_reads}_aln-to_{hifi_type}.{ont_type}.na.{chrom}.reads.fasta.gz',
        assembly = 'output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.na.{chrom}.fasta',
    output:
        chk = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.chk',
        bed = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}/{sample}_kmers_dist_diff.bed'
    log:
        'log/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.log'
    benchmark:
        'rsrc/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm.rsrc'
    wildcard_constraints:
        chrom = '(chrY|chrX)'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap_env"]}'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        walltime = lambda wildcards, attempt: f'{11*attempt}:59:00'
    params:
        preset = lambda wildcards: {'HIFIRW': 'hifi', 'ONTUL': 'ont'}[wildcards.other_reads],
        outdir = lambda wildcards, output: output.chk.rsplit('.', 2)[0]
    shell:
        'python repos/VerityMap/veritymap/main.py --no-reuse --reads {input.reads} -t {threads} '
        '-d {params.preset} -l {wildcards.sample} -o {params.outdir} {input.assembly} &> {log}'
            ' && '
        'touch {output.chk}'


rule normalize_veritymap_bed_file:
    """
    Normalize VerityMap BED file for simple
    merge with SNV/HET errors from variant
    callers - just format change...
    """
    input:
        bed = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}/{sample}_kmers_dist_diff.bed'
    output:
        tsv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.{other_reads}.vm-errors.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'est_size', 'num_reads'])
        df.drop('num_reads', axis=1, inplace=True)

        df['sample'] = wildcards.sample
        df['reads'] = wildcards.other_reads
        df['source'] = 'VerityMap'

        df.to_csv(output.tsv, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule merge_error_annotations:
    input:
        snv_dv = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.dv-errors.tsv',
        snv_pr = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.ONTUL.pr-errors.tsv',
        hifi_vm = 'output/eval/merged_errors/norm_tables/{sample}.{hifi_type}.{ont_type}.na.{chrom}.HIFIRW.vm-errors.tsv',
        quast = 'output/eval/assm_stats/SAMPLES.{hifi_type}.{ont_type}.na.{chrom}.quast-report.tsv'
    output:
        tsv = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.errors.tsv',
        stats = 'output/eval/merged_errors/{sample}.{hifi_type}.{ont_type}.na.{chrom}.error-stats.tsv',
    wildcard_constraints:
        chrom = '(chrX|chrY)'
    run:
        import pandas as pd

        assm_stats = pd.read_csv(input.quast, sep='\t', header=0)
        assm_length = assm_stats.at[assm_stats['sample'] == wildcards.sample, 'assembly_length_bp']

        all_errors = []
        error_stats = dict()
        for err_file in [input.snv_dv, input.snv_pr, input.hifi_vm]:
            df = pd.read_csv(err_file, sep='\t', header=0)
            all_errors.append(df)
        all_errors = pd.concat(all_errors, axis=0, ignore_index=False)
        all_errors.sort_values(['chrom', 'start', 'end'], ascending=True, inplace=True)

        error_stats['assembly_length_bp'] = assm_length
        error_stats['num_errors'] = all_errors.shape[0]
        error_stats['num_SNV_errors'] = all_errors.loc[all_errors['est_size'] == 1, :].shape[0]
        error_stats['error_bp'] = all_errors['est_size'].abs().sum()
        error_stats['error_bp_per_kbp'] = round(error_stats['bp_errors'] / assm_length * 1000, 5)
        error_stats['error_SNV_per_kbp'] = round(error_stats['num_SNV_errors'] / assm_length * 1000, 5)

        all_errors.to_csv(output.tsv, sep='\t', header=True, index=False)

        error_stats = pd.DataFrame(error_stats)
        error_stats.to_csv(output.stats, sep='\t', header=True, index=False)
    # END OF RUN BLOCK


###########################
# BELOW: DEPRECATED
# VerityMap cannot process
# a diploid genome
###########################


rule merge_read_sets:
    """
    VerityMap only supports a single
    input file for the reads...
    """
    input:
        reads = lambda wildcards: SAMPLE_DATA[wildcards.sample][wildcards.read_type]
    output:
        reads = temp('temp/merged_reads/{sample}_{read_type}.fa.gz')
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        walltime = lambda wildcards, attempt: f'{23*attempt}:59:00'
    shell:
        'seqtk seq -A -C {input.reads} | pigz -p {threads} --best > {output.reads}'


rule run_wg_veritymap:
    """
    NB: important to execute the main.py at the right
    location due to the setup horror explained above.
    Since main is messing with sys.path, VerityMap
    can be properly executed irrespective of the install/setup
    into "." locations by pip.
    """
    input:
        build_ok = 'repos/veritymap.build.ok',
        reads = 'temp/merged_reads/{sample}_{read_type}.fa.gz',
        assembly = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta',
    output:
        chk = 'output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.chk'
    log:
        'log/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.log'
    benchmark:
        'rsrc/output/eval/assm_errors/{sample}.{hifi_type}.{ont_type}.na.wg.{read_type}.vm.rsrc'
    singularity:
        f'{config["container_store"]}/{config["container"]["veritymap_env"]}'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 65536 * attempt,
        walltime = lambda wildcards, attempt: f'{23*attempt}:59:00'
    params:
        preset = lambda wildcards: {'HIFIRW': 'hifi', 'ONTUL': 'ont'}[wildcards.read_type],
        outdir = lambda wildcards, output: output.chk.rsplit('.', 2)[0]
    shell:
        'python repos/VerityMap/veritymap/main.py --no-reuse --reads {input.reads} -t {threads} '
        '-d {params.preset} -l {wildcards.sample} -o {params.outdir} {input.assembly} &> {log}'
            ' && '
        'touch {output.chk}'
