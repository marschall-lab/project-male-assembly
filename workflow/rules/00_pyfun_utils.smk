import pathlib
import json
import subprocess


def find_script_path(script_name):

    current_root = workflow.basedir  # location of main Snakefile
    assert current_root.endswith('workflow')

    script_folder = (current_root / pathlib.Path('scripts')).resolve(strict=True)

    all_py_scripts = script_folder.glob('**/*.py')
    this_py_script = [s for s in all_py_scripts if s.name == script_name]
    if len(this_py_script) != 1:
        raise RuntimeError(f'Cannot find requested script {script_name} starting at {script_folder}: {this_py_script}')
    return this_py_script[0]


def rsync(source, target):

    cmd = f'rsync --checksum {source} {target}'
    _ = subprocess.check_call(cmd, shell='/usr/bin/bash')
    return


def select_reference_genome(ref_name, fasta_index=False):

    available_references = config['reference_genomes']
    try:
        ref_genome = pathlib.Path(available_references[ref_name])
    except KeyError:
        raise ValueError(f'Requested reference "{ref_name}" is not available: {available_references}')

    if fasta_index:
        ref_genome = ref_genome.with_suffix('.fasta.fai')

    return ref_genome


def select_whole_genome_assembly(wildcards):

    assert hasattr(wildcards, 'sub_folder'), f'no sub_folder wildcard: {wildcards}'
    formatter = dict(wildcards)
    if wildcards.sub_folder == '00_raw':
        assm = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.na.wg/assembly.fasta'.format(**formatter)
    elif wildcards.sub_folder == '10_renamed':
        assm = 'output/hybrid/renamed/{sample}.{hifi_type}.{ont_type}.na.wg.fasta'.format(**formatter)
    else:
        raise ValueError(wildcards)
    return assm


def select_input_contig_alignment(wildcards):

    s_to_s = 'output/alignments/contigs-to-contigs/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{target}.paf.gz'
    s_to_ref = 'output/subset_wg/30_extract_ctgaln/{sample}.{hifi_type}.{ont_type}.na.chrY_aln-to_{target}.paf.gz'
    if sample.target in config['reference_genomes']:
        path = s_to_ref
    else:
        path = s_to_s
    path = path.format(**dict(wildcards))
    return path


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
                if major == 'bioconda':
                    # this is a release version w/o further
                    # "minor" qualification in the version
                    # string
                    major = version_info[-1]
                    minor = 'release'
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

    file_name = pl.Path(version_file).name
    assert file_name.endswith('.verkko.info')
    assembly_id = file_name.rsplit('.', 2)[0]

    major, minor = parse_verkko_version(version_file)
    assert major == config['verkko_major'].strip('"'), f'Verkko version error: {major} / {minor}'
    assert minor == config['verkko_minor'].strip('"'), f'Verkko version error: {major} / {minor}'
    subfolder = pl.Path(prefix, f'verkko_{major}_{minor}', suffix)
    full_path = share_path / subfolder
    full_path.mkdir(parents=True, exist_ok=True)
    return full_path, assembly_id


def get_fasta_seq_length(fasta_file):

    seq_lengths = dict()
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if not line.strip():
                continue
            if line.startswith('>'):
                seq_name = line.strip()[1:]
            else:
                seq_length = len(line.strip())
                seq_lengths[seq_name] = seq_length
    return seq_lengths
