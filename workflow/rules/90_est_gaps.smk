
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


rule generate_chry_graph_e2e_only:
    input:
        ref = "references_derived/T2T_chrY.fasta",
        assm = expand(
            "output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{{mapq}}.chrY.fasta",
            sample=E2E_SAMPLES,
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"]
        )
    output:
        ref_graph = "output/eval/par1_var/graphs/T2TY.2samples.{mapq}.chrY.gfa"
    log:
        "log/output/eval/par1_var/graphs/T2TY.2samples.{mapq}.chrY.gfa"
    conda:
        "../envs/graphtools.yaml"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 + 16384 * attempt,
        walltime = lambda wildcards, attempt: f'{11 ** attempt}:59:00'
    shell:
        "minigraph -t{threads} -cxggs {input.ref} {input.assm} "
        " > {output} 2> {log}"


localrules: create_graph_coloring
rule create_graph_coloring:
    input:
        ref_graph = "output/eval/par1_var/graphs/T2TY.{num_samples}samples.{mapq}.chrY.gfa"
    output:
        table = "output/eval/par1_var/graphs/T2TY.{num_samples}samples.{mapq}.chrY.annotations.csv"
    run:
        import pandas as pd
        # TODO: extract that from annotation sheet
        sample_infos = {
            "T2T": ("#808080", "J1a2a1a2c1a1", "ASK", "EUR"),
            "HG00358": ("#036C77", "N1a1a1a1a2a1a1a1a1a", "FIN", "EUR"),
            "HG02666": ("#D8A105", "A1a", "GWD", "AFR")
        }
        node_infos = []
        with open(input.ref_graph, "r") as gfa:
            for line in gfa:
                if not line.startswith("S"):
                    break
                columns = line.strip().split()
                node_id = columns[1]
                sample_id = columns[4]
                assert sample_id.startswith("SN")
                sample_id = sample.split(":")[-1]
                if sample_id == "chrY":
                    sample_id = "T2T"
                else:
                    sample_id = sample_id.split(".")[-1]
                    sample_id = sample_id.replace("HC", "HG")
                node_length = len(columns[2])
                if node_length == 1:
                    var_class = "SNV"
                elif node_length < 50:
                    var_class = "INDEL"
                else:
                    var_class = "SV"
                color, hapgroup, pop, spop = sample_infos[sample_id]
                node_infos.append((node_id, color, sample_id, hapgroup, pop, spop, var_class))
        df = pd.DataFrame.from_records(
            node_infos, columns=[
                "node", "color", "sample", "haplogroup", "population",
                "super_pop", "variant_type"
            ]
        )
        df.sort_values("node", inplace=True)
        df.to_csv(output.table, header=True, index=False)
    # END OF RUN BLOCK


rule run_all_graph_builds:
    input:
        gfa = expand(
            "output/eval/par1_var/graphs/T2TY.{sample_num}samples.{mapq}.chrY.gfa",
            mapq=["na"],
            sample_num=[10, 2]
        ),
