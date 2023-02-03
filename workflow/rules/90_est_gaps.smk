
rule merge_graph_infos:
    input:
        graph = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg/assembly.homopolymer-compressed.gfa',
        hifi_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg/assembly.hifi-coverage.csv',
        ont_cov = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg/assembly.ont-coverage.csv',
        y_names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.names.txt',
        x_names = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.chrX.names.txt',
    output:
        unitigs = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.unitigs.tsv",
        coloring = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.coloring.tsv",
        components = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.conncomp.tsv",
    conda:
        "../envs/pyscript.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
    params:
        script = find_script_path("merge_graph_infos.py"),
        layout = lambda wildcards, input: pl.Path(input.graph).parent.joinpath(
            "6-layoutContigs", "unitig-popped.layout.scfmap" 
        )
    shell:
        "{params.script} --graph {input.graph} --scf-layout {params.layout} "
        "--contig-names {input.y_names} {input.x_names} "
        "--hifi-node-cov {input.hifi_cov} --ont-node-cov {input.ont_cov} "
        "--out-unitigs {output.unitigs} --out-coloring {output.coloring} "
        "--out-conncomp {output.components}"


rule run_rukki:
    """ Stable Rukki was added to later releases
    of Verkko, hence this simply uses the latest
    Verkko release (v1.2)
    """
    input:
        graph = 'output/hybrid/verkko/{sample}.{hifi_type}.{ont_type}.{mapq}.wg/assembly.homopolymer-compressed.gfa',
        coloring = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.coloring.tsv",
    output:
        paths = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.rukki_paths.tsv",
        assign = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.rukki_assignments.tsv",
    singularity:
        f"{config['container_store']}/verkko_v1.2-release.sif"
    resources:
        mem_mb = lambda wildcards, attempt: 4096 + 2048 * attempt
    shell:
        "/repo/verkko/lib/verkko/bin/rukki trio "
        "-g {input.graph} -m {input.coloring} "
        "-p {output.paths} --final-assign {output.assign} "
        "--try-fill-bubbles --min-gap-size 1 --default-gap-size 1000"


rule add_rukki_paths:
    input:
        unitigs = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.unitigs.tsv",
        paths = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.rukki_paths.tsv",
    output:
        unitigs = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.unitig-paths.tsv",
        subset = "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.paths-chrY.tsv",
    conda:
        "../envs/pyscript.yaml"
    params:
        script = find_script_path("add_rukki_paths.py")
    shell:
        "{params.script} -u {input.unitigs} -p {input.paths} "
        "--rukki-default-gap 1000 --summary-subset chrY "
        "--out-table {output.unitigs} --out-summary {output.subset}"


rule run_all_gap_estimates:
    input:
        subsets = expand(
            "output/eval/assm_gaps/proc/{sample}.{hifi_type}.{ont_type}.{mapq}.paths-chrY.tsv",
            sample=COMPLETE_SAMPLES,
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            mapq=["na"]
        )


rule generate_chry_graph:
    input:
        ref = "references_derived/T2T_chrY.fasta",
        assm = expand(
            "output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{{mapq}}.chrY.fasta",
            sample=PAR1_SAMPLES,
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"]
        )
    output:
        ref_graph = "output/eval/par1_var/graphs/T2TY.10samples.{mapq}.chrY.gfa"
    log:
        "log/output/eval/par1_var/graphs/T2TY.10samples.{mapq}.chrY.gfa"
    conda:
        "../envs/graphtools.yaml"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 + 16384 * attempt,
        walltime = lambda wildcards, attempt: f'{11 ** attempt}:59:00'
    shell:
        "minigraph -t{threads} -cxggs {input.ref} {input.assm} "
        " > {output} 2> {log}"


rule run_all_graph_builds:
    input:
        gfa = expand(
            "output/eval/par1_var/graphs/T2TY.10samples.{mapq}.chrY.gfa",
            mapq=["na"]
        )
