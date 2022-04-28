
import pathlib
import pandas


def read_sample_table():

    debug = bool(config.get('debug', False))

    if not 'samples' in config:
        raise ValueError('No key "samples" in config - should be a path to the sample table (TSV')
    sample_table_file = pl.Path(config['samples']).resolve(strict=True)
    samples = pd.read_csv(
        sample_table_file,
        sep='\t',
        header=0
    )
    mandatory_columns = ['sample', 'hifi', 'ont']
    assert all(c in samples.columns for c in mandatory_columns), \
        f'Sample table does not contain all mandatory columns: {mandatory_columns}'

    selector = []
    for row in samples.itertuples(index=False):
        try:
            _ = pathlib.Path(row.hifi).resolve(strict=True)
        except FileNotFoundError:
            if debug:
                sys.stderr.write(f'\nError skip: no (valid) HiFi path for sample {sample}\n')
            selector.append(False)
            continue
        try:
            _ = pathlib.Path(row.ont).resolve(strict=True)
        except FileNotFoundError:
            if debug:
                sys.stderr.write(f'\nError skip: no (valid) ONT path set for sample {sample}\n')
            selector.append(False)
            continue
        selector.append(True)

    samples = samples.loc[selector, :].copy()
    return samples


def collect_sample_data(samples):

    sample_data = collections.defaultdict(dict)

    # NB: by construction, samples contains
    # only complete samples (HiFi + ONT available)
    for row in samples.itertuples(index=False):
        hifi_data = get_hifi_data(row.hifi)
        sample_data[row.sample]['HIFIRW'] = hifi_data
        sample_data[row.sample][('HIFIRW', 'num')] = len(hifi_data)
        ont_data = get_ont_data(row.ont)
        sample_data[row.sample]['ONTUL'] = ont_data
        sample_data[row.sample][('ONTUL', 'num')] = len(ont_data)
    return sample_data


def get_hifi_data(hifi_path):

    suffix = '.fastq.gz'
    hifi_files = get_data(hifi_path, suffix)
    return hifi_files


def get_ont_data(ont_path):

    suffix = 'guppy-5.0.11-sup-prom_fastq_pass.fastq.gz'
    ont_files = get_data(ont_path, suffix)
    return ont_files


def get_data(top_path, suffix):

    search_path = pathlib.Path(top_path).resolve(strict=True)
    data_files = sorted(search_path.glob(f'**/{suffix}'))
    assert len(data_files) > 0, f'No data found: {suffix} / {top_path}'
    return data_files


# BELOW: OLD DATA COLLECTION FUNCTIONS
# DELETE WHEN FREEZE 1 IS AVAILABLE


# def read_1000gp_sample_table(tsv_path):
#     """
#     This function reads a sample table
#     containing all 1000 GP samples, and
#     extracts all males.
#     """
#     import pandas as pd
#     import collections as col

#     df = pd.read_csv(tsv_path, sep='\t', header=0)
#     sample_infos = col.defaultdict(dict)

#     for row in df.itertuples(index=False):
#         sex = row.sex[0].upper()
#         assert sex in ['M', 'F']
#         long_sample = f'{row.super_population}-{row.population}-{row.family_id}-{sex}_{row.sample}'
#         sample_infos[row.sample]['long_id'] = long_sample
#         sample_infos[row.sample]['sex'] = sex

#     # add several special samples for merging reads
#     sample_infos['NA193N7']['long_id'] = 'AFR-LWK-DUO-M_NA193N7'
#     sample_infos['NA193N7']['sex'] = 'M'
#     sample_infos['NA193N7']['merge'] = ['NA19317', 'NA19347']
#     sample_infos['NA193N7']['is_regular'] = False

#     sample_infos['NA193NN']['long_id'] = 'AFR-LWK-TRIO-M_NA193NN'
#     sample_infos['NA193NN']['sex'] = 'M'
#     sample_infos['NA193NN']['merge'] = ['NA19317', 'NA19347', 'NA19384']
#     sample_infos['NA193NN']['is_regular'] = False

#     sample_infos['AFR4MIX']['long_id'] = 'AFR-MIX-QUART-M_AFR4MIX'
#     sample_infos['AFR4MIX']['sex'] = 'M'
#     sample_infos['AFR4MIX']['merge'] = ['NA19317', 'NA19347', 'NA19384', 'HG02666']
#     sample_infos['AFR4MIX']['is_regular'] = False

#     return sample_infos


# def normalize_sample_name(sample_infos, test_sample):

#     try:
#         if test_sample in sample_infos:
#             alt_sample = test_sample
#         elif 'NA' in test_sample:
#             alt_sample = 'GM' + test_sample[2:]
#         elif 'GM' in test_sample:
#             alt_sample = 'NA' + test_sample[2:]
#         else:
#             alt_sample = 'undetermined'
#         _ = sample_infos[alt_sample]
#     except KeyError:
#         raise KeyError(f'Sample / ALT sample name does not exist: {test_sample} / {alt_sample}')
#     return alt_sample


# def add_hifirw_readsets(sample_infos, hifirw_path):

#     suffix = '.fastq.gz'
#     flag_files = pathlib.Path(hifirw_path).glob('**/*.final')
#     hifirw_samples = []

#     for flag_file in flag_files:
#         assert flag_file.name == 'status.final'
#         sample = flag_file.parent.name
#         # raw / data not preprocessed via QC pipeline
#         # may not yet have a normalized sample name
#         sample = normalize_sample_name(sample_infos, sample)

#         fastq_files = sorted(flag_file.parent.glob(f'*{suffix}'))
#         if not fastq_files:
#             raise ValueError(f'No HiFi (raw) FASTQ (split) files for sample / path: {sample} / {flag_file.parent}')
#         sample_infos[sample]['HIFIRW'] = fastq_files
#         hifirw_samples.append(sample)

#     return sample_infos, sorted(hifirw_samples)


# def add_hifiec_readsets(sample_infos, hifiec_path):

#     fasta_files = pathlib.Path(hifiec_path).glob('*.ec-reads.fasta.gz')
#     hifiec_samples = []

#     for fasta_file in fasta_files:
#         file_size_bytes = fasta_file.stat().st_size
#         if file_size_bytes < 15 * 1024**3:
#             # assume file is still being dumped to disk
#             continue
#         file_sample = fasta_file.name.split('_')[0]
#         assert file_sample in sample_infos
#         sample_infos[file_sample]['HIFIEC'] = fasta_file
#         hifiec_samples.append(file_sample)
#     return sample_infos, sorted(hifiec_samples)


# def add_hifiaf_readsets(sample_infos, hifiaf_path):

#     fastq_files = pathlib.Path(hifiaf_path).glob('*_1000.fastq.gz')
#     hifiaf_samples = []

#     for fastq_file in fastq_files:
#         file_size_bytes = fastq_file.stat().st_size
#         if file_size_bytes < 15 * 1024**3:
#             # assume file is still being dumped to disk
#             continue
#         file_sample = fastq_file.name.split('_')[0]
#         assert file_sample in sample_infos
#         sample_infos[file_sample]['HIFIAF'] = fastq_file
#         hifiaf_samples.append(file_sample)
#     return sample_infos, sorted(hifiaf_samples)


# def add_assembly_graphs(sample_infos, assembly_path):
#     import collections as col

#     graph_files = pathlib.Path(assembly_path).glob('**/*.gfa')
#     assembled_samples = set()
#     graph_count = col.Counter()

#     for graph_file in graph_files:
#         if 'noseq' in graph_file.name:
#             continue
#         if 'p_utg' in graph_file.name:
#             continue
#         file_sample = graph_file.name.split('_')[0]
#         assert file_sample in sample_infos
#         if 'a_ctg' in graph_file.name:
#             sample_infos[file_sample]['TIGALT'] = graph_file
#         elif 'p_ctg' in graph_file.name:
#             sample_infos[file_sample]['TIGPRI'] = graph_file
#         elif 'r_utg' in graph_file.name:
#             sample_infos[file_sample]['TIGRAW'] = graph_file
#         else:
#             raise ValueError(f'Unexpected file: {graph_file.name}')
#         assembled_samples.add(file_sample)
#         graph_count[file_sample] += 1
#     for n, c in graph_count.most_common():
#         if c != 3:
#             raise ValueError(f'Missing assembly graph type for sample {n}')
#     return sample_infos, sorted(assembled_samples)


# def add_ontul_readsets_merged(sample_infos, ontul_path):

#     suffix = 'GPYv5011SUP.fasta.gz'
#     fasta_files = pathlib.Path(ontul_path).glob(f'*{suffix}')
#     ontul_samples = []

#     for fasta_file in fasta_files:
#         file_size_bytes = fasta_file.stat().st_size
#         if file_size_bytes < 30 * 1024**3:
#             # assume file is still incomplete
#             continue
#         file_sample = fasta_file.name.split('_')[0]
#         assert file_sample in sample_infos
#         sample_infos[file_sample]['ONTUL'] = fasta_file
#         ontul_samples.append(file_sample)

#     return sample_infos, sorted(ontul_samples)


# def add_ontul_readsets_split(sample_infos, ontul_path):

#     suffix = 'guppy-5.0.11-sup-prom_fastq_pass.fastq.gz'
#     flag_files = pathlib.Path(ontul_path).glob('**/*.final')
#     ontul_samples = []

#     for flag_file in flag_files:
#         assert flag_file.name == 'status.final'
#         sample = flag_file.parent.name
#         # raw / data not preprocessed via QC pipeline
#         # may not yet have a normalized sample name
#         sample = normalize_sample_name(sample_infos, sample)

#         fastq_files = sorted(flag_file.parent.glob(f'*{suffix}'))
#         if not fastq_files:
#             raise ValueError(f'No ONT-UL FASTQ (split) files for sample / path: {sample} / {flag_file.parent}')
#         sample_infos[sample]['ONTUL'] = fastq_files
#         ontul_samples.append(sample)

#     return sample_infos, sorted(ontul_samples)


# def add_ontec_readsets(sample_infos, ontec_path):

#     suffix = 'ONTEC_B9A766C4.fasta.gz'
#     # hash key = combination of read sets used to generate
#     # this ONTEC read set
#     fasta_files = pathlib.Path(ontec_path).glob(f'*{suffix}')
#     ontec_samples = []

#     for fasta_file in fasta_files:
#         # since the ONTEC FASTA files are just symlinked
#         # into the input folder, no size check necessary
#         file_sample = fasta_file.name.split('_')[0]
#         assert file_sample in sample_infos
#         sample_infos[file_sample]['ONTEC'] = fasta_file
#         ontec_samples.append(file_sample)

#     return sample_infos, sorted(ontec_samples)


# PATH_SAMPLE_TABLE = config['path_sample_table']
# PATH_HIFIRW_READS = config['path_hifirw_reads']
# PATH_HIFIEC_READS = config['path_hifiec_reads']
# PATH_HIFIAF_READS = config['path_hifiaf_reads']
# PATH_ASSEMBLY_GRAPHS = config['path_assembly_graphs']
# if config.get('use_preprocessed_readsets', False):
#     PATH_ONTUL_READS = config['path_ontul_reads_merged']
# else:
#     PATH_ONTUL_READS = config['path_ontul_reads_split']
# PATH_ONTEC_READS = config['path_ontec_reads']


# def init_samples_and_data():

#     location_smk_file = pathlib.Path(workflow.basedir)
#     location_sample_table = (location_smk_file / pathlib.Path(PATH_SAMPLE_TABLE)).resolve()
    
#     sample_infos = read_sample_table(location_sample_table)
#     hifirw_samples = []
#     for hifi_path in PATH_HIFIRW_READS:
#         sample_infos, path_hifirw_samples = add_hifirw_readsets(sample_infos, hifi_path)
#         hifirw_samples.extend(path_hifirw_samples)
#     sample_infos, hifiec_samples = add_hifiec_readsets(sample_infos, PATH_HIFIEC_READS)
#     sample_infos, hifiaf_samples = add_hifiaf_readsets(sample_infos, PATH_HIFIAF_READS)
#     sample_infos, assembled_samples = add_assembly_graphs(sample_infos, PATH_ASSEMBLY_GRAPHS)
#     if config.get('use_preprocessed_readsets', False):
#         sample_infos, ontul_samples = add_ontul_readsets_merged(sample_infos, PATH_ONTUL_READS)
#     else:
#         ontul_samples = []
#         for ont_path in PATH_ONTUL_READS:
#             sample_infos, path_ont_samples = add_ontul_readsets_split(sample_infos, ont_path)
#             ontul_samples.extend(path_ont_samples)
#     sample_infos, ontec_samples = add_ontec_readsets(sample_infos, PATH_ONTEC_READS)

#     return sample_infos, hifirw_samples, hifiec_samples, hifiaf_samples, ontul_samples, ontec_samples, assembled_samples


# SAMPLE_INFOS, HIFIRW_SAMPLES, HIFIEC_SAMPLES, HIFIAF_SAMPLES, ONTUL_SAMPLES, ONTEC_SAMPLES, ASSEMBLED_SAMPLES = init_samples_and_data()

# CONSTRAINT_REGULAR_SAMPLES = '(' + '|'.join(sorted(k for k in SAMPLE_INFOS.keys() if SAMPLE_DATA[k].get('is_regular', True))) + ')'
# CONSTRAINT_ALL_SAMPLES = '(' + '|'.join(sorted(SAMPLE_INFOS.keys())) + ')'
