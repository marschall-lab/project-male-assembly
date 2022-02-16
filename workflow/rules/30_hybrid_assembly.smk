import pathlib

localrules: run_verkko_test_data, run_verkko_whole_genome_assembly

rule run_verkko_test_data:
    input:
        hifi = '/gpfs/project/projects/medbioinf/data/share/globus/dev_debug/ga_mbg/verkko_test-data/hifi.fastq.gz',
        ont = '/gpfs/project/projects/medbioinf/data/share/globus/dev_debug/ga_mbg/verkko_test-data/ont.fastq.gz',
        verkko_install = '/gpfs/project/ebertp/data/repository/verkko'  # install directory w/ code changes
    output:
        assembly = 'output/temp/test_verkko/assembly.fasta'
    log:
        'log/output/temp/test_verkko/run.log'
    conda:
        '../envs/verkko_test.yaml'
    params:
        verkko_bin = lambda wildcards, input: pathlib.Path(input.verkko_install, 'bin').resolve(strict=True),
        verkko_lib_bin = lambda wildcards, input: pathlib.Path(input.verkko_install, 'lib/verkko/bin').resolve(strict=True),
        workdir = lambda wildcards, output: pathlib.Path(output.assembly).parent
    shell:
        'module load gcc/10.2.0 ; '
        'VERKKO={input.verkko_install} '
        'PATH={params.verkko_bin}:{params.verkko_lib_bin}:$PATH '
        'verkko -d {params.workdir} --hifi {input.hifi} --nano {input.ont} '
        '--python `which python` --mbg `which MBG` --graphaligner `which GraphAligner` &> {log}'


rule run_verkko_targeted_assembly:
    """
    Note that the below call explicitly deletes the working
    directory before starting Verkko as Verkko's own
    Snakemake workflow cannot necessarily detect corrupt
    leftovers from previous runs if these were caused by
    bugs in one of the scripts.
    In principle, Verkko should be capable to continue
    the previous run, but I leave that for a post-beta
    release version for now.
    """
    input:
        ont = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{ont_type}.{chrom}-reads.{mapq}.fasta.gz',
        hifi = 'output/read_subsets/{chrom}/{sample_info}_{sample}_{hifi_type}.{chrom}-reads.{mapq}.fasta.gz'
    output:
        assembly = multiext(
            'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}/assembly',
            '.fasta', '.gfa', '.hifi-coverage.csv', '.layout', '.ont-coverage.csv'
        ),
        version = 'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verrko.info'
    log:
        'log/output/hybrid/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.log'
    benchmark:
        'rsrc/output/hybrid/{sample_info}_{sample}.{hifi_type}.{ont_type}.{mapq}.{chrom}.verkko.rsrc'
    wildcard_constraints:
        chrom = '(chrY|chrX|chrXY)'
    # conda:
    #     '../envs/verkko.yaml'
    singularity:
        'verkko.sif'
    threads: config['num_cpu_high']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * (12288 if wildcards.hifi_type in ['HIFIEC', 'ONTEC', 'OHEC'] else 24576),
        walltime = lambda wildcards, attempt: f'{5 ** attempt}:59:00'
    params:
        run_correction = lambda wildcards: '--no-correction' if wildcards.hifi_type in ['HIFIEC', 'ONTEC', 'OHEC'] else '',
        high_maxk = lambda wildcards: '--max-k 100000' if wildcards.hifi_type in ['ONTEC', 'OHEC'] else '',
        work_dir = lambda wildcards, output: output.assembly[0].rsplit('/', 1)[0]
    shell:
        'rm -rf {params.work_dir} && '
        'cat /.singularity.d/labels.json > {output.version} && '
        '/repos/verkko/bin/verkko --local --hifi {input.hifi} --nano {input.ont} -d {params.work_dir} '
            '{params.run_correction} {params.high_maxk} --threads {threads} &> {log} '


rule run_verkko_whole_genome_assembly:
    """
    Whole-genome runs are only feasible on a distributed cluster infrastructure,
    which requires some code changes to Verkko for the HILBERT cluster. This mimicks
    the rule above for the test data.

    Note about Verkko / Canu scripts:
    - compiled w/ gcc/10.2.0

    Note about MBG dev version:
    - compiled w/ gcc/10.2.0 and zlib/1.2.11
    - export LIBRARY_PATH=$LD_LIBRARY_PATH

    """
    input:
        hifi = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.hifi_type],
        ont = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type],
        verkko_install = '/gpfs/project/ebertp/data/repository/verkko',
        mbg_install = '/gpfs/project/ebertp/data/repository/MBG' 
    output:
        assembly = multiext(
            'output/hybrid/verkko/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg/assembly',
            '.fasta', '.gfa', '.hifi-coverage.csv', '.layout', '.ont-coverage.csv'
        ),
    log:
        'log/output/hybrid/{sample_info}_{sample}.{hifi_type}.{ont_type}.na.wg.verkko.log'
    conda:
        '../envs/verkko_test.yaml'
    params:
        verkko_bin = lambda wildcards, input: pathlib.Path(input.verkko_install, 'bin').resolve(strict=True),
        verkko_lib_bin = lambda wildcards, input: pathlib.Path(input.verkko_install, 'lib/verkko/bin').resolve(strict=True),
        mbg_bin = lambda wildcards, input: pathlib.Path(input.mbg_install, 'bin').resolve(strict=True),
        workdir = lambda wildcards, output: pathlib.Path(output.assembly).parent
    shell:
        'module load gcc/10.2.0 ; '
        'module load zlib/1.2.11 ; '
        'VERKKO={input.verkko_install} '
        'PATH={params.verkko_bin}:{params.verkko_lib_bin}:$PATH '
        'verkko -d {params.workdir} --hifi {input.hifi} --nano {input.ont} '
        '--python `which python` --mbg {params.mbg_bin}/MBG --graphaligner `which GraphAligner` &> {log}'
