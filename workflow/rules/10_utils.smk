import hashlib

localrules: dump_mbg_param_info


MBG_PARAMS = compute_lut_mbg_params()
LJA_PARAMS = compute_lut_lja_params()
HIFIASM_TIGS = config['hifiasm_tig_names']


rule dump_mbg_param_info:
    output:
        expand(
            'assembler_params/MBG_{param_info}.info',
            param_info=['{}_k{}-w{}-r{}'.format(phash, *pvalues) for phash, pvalues in MBG_PARAMS.items()]
        )
    priority: 100
    run:
        file_template = 'assembler_params/MBG_{}.info'
        for phash, pvalues in MBG_PARAMS.items():
            param_info = '{}_k{}-w{}-r{}'.format(phash, *pvalues)
            with open(file_template.format(param_info), 'w') as dump:
                pass


rule index_bam_alignment:
    input:
        '{filepath}.bam'
    output:
        '{filepath}.bam.bai'
    conda:
        '../envs/biotools.yaml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        walltime = lambda wildcards, attempt: f'{attempt ** 3:02}:59:00'
    shell:
        'samtools index -b -@ {threads} {input}'
