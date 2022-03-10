import pathlib
import json


def find_script_path(script_name, subfolder=''):
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-male-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


def select_alignment_cache_file(wildcards):

    input_qc_reference = config['input_qc_reference']

    cache_dir = pl.Path(config['path_alignment_cache'])
    glob_expression = f'{wildcards.sample}*{wildcards.read_type}*{input_qc_reference}.cov.cache.h5'
    cache_files = list(cache_dir.glob(glob_expression))
    if not cache_files:
        raise ValueError(
            'Missing input - no cached whole-genome read coverage file for expression '
            f'{glob_expression} /// {str(wildcards)}'
        )
    if len(cache_files) > 1:
        raise ValueError(
            'Ambiguous input - more than one cached whole-genome read coverage file for expression '
            f'{glob_expression} /// {str(wildcards)}: {cache_files}'
        )
    return cache_files[0]


def select_afr_mix_subsets(wildcards):

    assert 'mq0' in wildcards.mapq
    merge_samples = SAMPLE_INFOS[wildcards.sample]['merge']
    template = 'output/read_subsets/{chrom}/{sample_long}_{read_type}.{chrom}-reads.{mapq}.fasta.gz'
    selected_readsets = []
    for sample in merge_samples:
        sample_long = SAMPLE_INFOS[sample]['long_id']
        formatter = {
            'sample_long': sample_long,
            'read_type': wildcards.read_type,
            'mapq': wildcards.mapq,
            'chrom': wildcards.chrom
        }
        selected_readsets.append(template.format(**formatter))
    return sorted(selected_readsets)


def select_hybrid_assembly_input_graph(wildcards):
    assemblers = {
        'HAS': 'hifiasm',
        'MBG': 'mbg',
        'LJA': 'lja'
    }
    graphs_paths = {
        'hifiasm': 'output/target_assembly/{chrom}/hifiasm/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.{tigs}.gfa',
        'mbg': 'output/target_assembly/{chrom}/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.gfa',
        'lja': 'no_uncompressed_graph_available'
    }
    try:
        # for hifiasm whole-genome graphs
        gfa = SAMPLE_INFOS[wildcards.sample][wildcards.tigs]
    except KeyError:
        if 'OHEC' in wildcards.tigs:
            read_type = 'OHEC'
        elif 'OEC' in wildcards.tigs:
            read_type = 'ONTEC'
        elif 'EC' in wildcards.tigs:
            read_type = 'HIFIEC'
        elif 'AF' in wildcards.tigs:
            read_type = 'HIFIAF'
        else:
            raise
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise
        assembler = assemblers[wildcards.tigs[:3]]
        if assembler == 'lja':
            raise ValueError(f'Not supported at the moment - assembler LJA / {str(wildcards)}')
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'read_type': read_type,
            'mapq': mapq,
            'chrom': chrom
        }
        if assembler == 'hifiasm':
            formatter['tigs'] = get_hifiasm_tigs(wildcards.tigs)
        elif assembler == 'mbg':
            k, w, r = get_mbg_param(wildcards.tigs)
            formatter['kmer'] = k
            formatter['window'] = w
            formatter['resolvek'] = r
        else:
            raise
        gfa = graphs_paths[assembler].format(**formatter)
    return gfa


def select_hybrid_assm_ont_reads(wildcards):

    if 'MQ' in wildcards.tigs:
        assert wildcards.ont_type == 'ONTUL'
        assert 'MQ0' in wildcards.tigs
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise ValueError(f'Unknown chrom: {str(wildcards)}')
        template = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTUL.{chrom}-reads.{mapq}.fasta.gz'
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'mapq': mapq,
            'chrom': chrom
        }
        ont_reads = template.format(**formatter)
    else:
        ont_reads = str(SAMPLE_INFOS[wildcards.sample][wildcards.ont_type])
    if 'gpfs' in ont_reads:
        ont_reads = str(pl.Path('/hilbert', ont_reads.strip('/')))
    return ont_reads


def select_hybrid_assm_hifi_reads(wildcards):

    if 'MQ' in wildcards.tigs:
        assert wildcards.ont_type == 'ONTUL'
        assert 'MQ0' in wildcards.tigs
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise ValueError(f'Unknown chrom: {str(wildcards)}')
        template = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTUL.{chrom}-reads.{mapq}.fasta.gz'
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'mapq': mapq,
            'chrom': chrom
        }
        ont_reads = template.format(**formatter)
    else:
        ont_reads = str(SAMPLE_INFOS[wildcards.sample][wildcards.ont_type])
    if 'gpfs' in ont_reads:
        ont_reads = str(pl.Path('/hilbert', ont_reads.strip('/')))
    return ont_reads


def get_hifiasm_tigs(tig_spec):

    try:
        result = HIFIASM_TIGS[tig_spec]
    except KeyError:
        try:
            result = HIFIASM_TIGS[tig_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown MBG param spec: {tig_spec} / {HIFIASM_TIGS}')
    return result


def compute_lut_mbg_params():

    init_kmer = config['mbg_init_kmer']
    window_sizes = config['mbg_window_size']
    resolve_kmer = config['mbg_resolve_kmer']

    if not len(init_kmer) == len(window_sizes) == len(resolve_kmer):
        raise ValueError('The same number of parameters has to be specified for MBG k/w/r')

    param_lut = dict()
    for k, w, r in zip(init_kmer, window_sizes, resolve_kmer):
        param_comb = f'k{k}-w{w}-r{r}'.encode('ascii')
        param_hash = hashlib.sha1(param_comb).hexdigest()
        key = param_hash.upper()[:6]
        if key in param_lut:
            raise ValueError(f"MBG param key collision: {param_comb.decode('ascii')} / {param_lut[key]}")
        param_lut[key] = k, w, r
    return param_lut


def get_mbg_param(param_spec, which=None):

    param_pos = {
        'k': 0,
        'w': 1,
        'r': 2,
        'kmer': 0,
        'window': 1,
        'resolve': 2
    }

    try:
        result = MBG_PARAMS[param_spec]
    except KeyError:
        try:
            result = MBG_PARAMS[param_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown MBG param spec: {param_spec} / {MBG_PARAMS}')
    if which is not None:
        assert which in param_pos, f'Cannot process "which": {which}'
        result = result[param_pos[which]]
    return result


def compute_lut_lja_params():

    small_kmer = config['lja_small_kmer']
    large_kmer = config['lja_large_kmer']

    if not len(small_kmer) == len(large_kmer):
        raise ValueError('The same number of parameters has to be specified for LJA k/K')

    param_lut = dict()
    for k, K in zip(small_kmer, large_kmer):
        param_comb = f'k{k}-K{K}'.encode('ascii')
        param_hash = hashlib.sha1(param_comb).hexdigest()
        key = param_hash.upper()[:6]
        if key in param_lut:
            raise ValueError(f"LJA param key collision: {param_comb.decode('ascii')} / {param_lut[key]}")
        param_lut[key] = k, K
    return param_lut


def get_lja_param(param_spec, which=None):

    param_pos = {
        'k': 0,
        'K': 1,
        'smallk': 0,
        'largek': 1
    }

    try:
        result = LJA_PARAMS[param_spec]
    except KeyError:
        try:
            result = LJA_PARAMS[param_spec.split('-')[1]]
        except KeyError:
            raise ValueError(f'Unknown LJA param spec: {param_spec} / {LJA_PARAMS}')
    if which is not None:
        assert which in param_pos, f'Cannot process "which": {which}'
        result = result[param_pos[which]]
    return result


def select_reference_genome(ref_name, fasta_index=False):

    available_references = config['reference_genomes']
    try:
        ref_genome = pathlib.Path(available_references[ref_name])
    except KeyError:
        raise ValueError(f'Requested reference "{ref_name}" is not available: {available_references}')

    if fasta_index:
        ref_genome = ref_genome.with_suffix('.fasta.fai')

    return ref_genome


def parse_verkko_version(file_path):

    major = None
    minor = None
    
    with open(file_path, 'r') as info_file:
        try:
            version = json.load(info_file)
            # custom build in Singularity Container
            # not suitable for whole-genome runs, thus
            # deprecated
            major = version['verkko_release']
            minor = version['verkko_commit']
        except json.JSONDecodeError:
            info_file.seek(0)
            version_info = info_file.read().strip().split()
            if 'bioconda' in version_info:
                major = version_info[-2]
                minor = version_info[-1]
            elif 'HEAD' in version_info:
                # git dev build, has no major version,
                # so use date timestamp
                major = 'dev-' + version_info[-3].replace('-', '')
                minor = version_info[-4][:7]
            else:
                raise ValueError(f'Unexpected Verkko version info: {version_info}')
    return major, minor


def determine_verkko_subfolder(version_file, prefix, suffix):
    """
    Create versioned Verkko subfolder for sharing results
    """
    share_path = pl.Path(config['path_root_share_working']).resolve(strict=True)

    file_name = pl.Path(versio_file).name
    assert file_name.endswith('.verkko.info')
    assembly_id = file_name.rsplit('.', 2)[0]

    major, minor = parse_verkko_version(versio_file)
    assert major == config['verkko_major'].strip('"'), f'Verkko version error: {major} / {minor}'
    assert minor == config['verkko_minor'].strip('"'), f'Verkko version error: {major} / {minor}'
    subfolder = pl.Path(prefix, f'verkko_{major}_{minor}', suffix)
    full_path = share_path / subfolder
    full_path.mkdir(parents=True, exist_ok=True)
    return full_path, assembly_id
