# Project: chromosome Y assembly

This project repository contains a Snakemake workflow to produce whole-genome Verkko
assemblies, and extract contigs that most likely represent the Y chromosome.
The workflow requires HiFi and ONT reads to be executed, plus Illumina short reads for
certain assembly evaluation tasks.

The input sample sheet is a simple tab-separated table listing sample ID and (file system)
location of read sets. The sample sheet needs to be loaded as follows:

```bash
$ snakemake --config samples=PATH_TO_SAMPLE_SHEET [...]
```

[Example sample sheet](data/samples.tsv)

## Tool setup

The entire workflow uses Conda environments wherever possible to deploy software dependencies.
A base environment containing Snakemake itself is defined in
[`workflow/envs/run_env.yaml`](workflow/envs/run_env.yaml).

For software in development/prototype stage (Verkko and VerityMap), adaptations to the local
infrastructure (Verkko) or building specific bugfix versions (VerityMap, see module
`workflow/envs/80_est_assm_errors.smk`) is required, with the former not being automatable.

## Plotting

The folder `notebooks/` contains Jupyter notebooks used to plot various summary statistics
of the generated assemblies. The notebooks contain a brief description documenting the necessary
input files (produced by the Snakemake workflow).


## Citation

In preparation
