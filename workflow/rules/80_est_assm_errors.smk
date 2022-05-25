
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

localrules: clone_veritymap_repo, build_veritymap, clean_veritymap_build


rule clone_veritymap_repo:
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
        'git checkout v2.1.1-alpha'
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


rule run_veritymap:
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
