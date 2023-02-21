
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
        assm = lambda wildcards: expand(
            "output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta",
            sample=GRAPH_SAMPLES[int(wildcards.num)],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            mapq=["na"]
        )
    output:
        ref_graph = "output/eval/par1_var/graphs/T2T.{num}samples.na.chrY.L{min_var_size}kbp.gfa"
    log:
        "log/output/eval/par1_var/graphs/T2T.{num}samples.na.chrY.L{min_var_size}kbp.minigraph.log"
    benchmark:
        "rsrc/output/eval/par1_var/graphs/T2T.{num}samples.na.chrY.L{min_var_size}kbp.minigraph.rsrc"
    conda:
        "../envs/graphtools.yaml"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 16384 + 16384 * attempt,
        walltime = lambda wildcards, attempt: f'{11 ** attempt}:59:00'
    params:
        min_var_size = lambda wildcards: int(wildcards.min_var_size) * 1000
    shell:
        "minigraph -t{threads} -L {params.min_var_size} -cxggs {input.ref} {input.assm} "
        " > {output} 2> {log}"


rule generate_chrxy_graph:
    input:
        ref_y = "references_derived/T2T_chrY.fasta",
        ref_x = "references_derived/T2T_chrX.fasta",
        assm_y = lambda wildcards: expand(
            "output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta",
            sample=GRAPH_SAMPLES[int(wildcards.num)],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            mapq=["na"]
        ),
        assm_x = lambda wildcards: expand(
            "output/subset_wg/20_extract_contigs/{sample}.{hifi_type}.{ont_type}.{mapq}.chrX.fasta",
            sample=GRAPH_SAMPLES[int(wildcards.num)],
            hifi_type=["HIFIRW"],
            ont_type=["ONTUL"],
            mapq=["na"]
        )
    output:
        ref_graph = "output/eval/par1_var/graphs/T2T.{num}samples.na.chrXY.L{min_var_size}kbp.gfa"
    log:
        "log/output/eval/par1_var/graphs/T2T.{num}samples.na.chrXY.L{min_var_size}kbp.minigraph.log"
    benchmark:
        "rsrc/output/eval/par1_var/graphs/T2T.{num}samples.na.chrXY.L{min_var_size}kbp.minigraph.rsrc"
    conda:
        "../envs/graphtools.yaml"
    threads: config['num_cpu_medium']
    resources:
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        walltime = lambda wildcards, attempt: f'{11 ** attempt}:59:00'
    shell:
        "minigraph -t{threads} -L {params.min_var_size} -cxggs "
        "{input.ref_y} {input.assm_y} {input.ref_x} {input.assm_x} "
        " > {output} 2> {log}"


localrules: prepare_graph_coloring
rule prepare_graph_coloring:
    input:
        color_table = config["colors"]
    output:
        dump = "references_derived/sample_annotations.json"
    run:
        import json
        import pandas as pd

        df = pd.read_csv(input.color_table, header=0, sep="\t")
        annotations = dict()
        for row in df.itertuples(index=False):
            if row.sample == "NA24385/HG002":
                sample = "T2T.Y"
                display_sample = "T2T.Y"
                display_color = "#D3D3D3"  # light grey
            else:
                sample = f"{row.workflow_sample}.Y"
                display_sample = f"{row.sample}.Y"
                display_color = row.superpop_color_hex
            infos = {
                "color": display_color,
                "haplogroup": row.Y_haplogroup,
                "population": row.population,
                "super_pop": row.super_population,
                "sample": display_sample
            }
            annotations[sample] = dict(infos)
            infos["haplogroup"] = "female"
            infos["sample"] = infos["sample"].replace(".Y", ".X")
            if infos["sample"] == "T2T.X":
                infos["color"] = "#5A5A5A"  # dark grey
            annotations[sample.replace(".Y", ".X")] = dict(infos)
        with open(output.dump, "w") as dump:
            json.dump(annotations, dump, ensure_ascii=True, indent=2)
    # END OF RUN BLOCK


localrules: create_graph_coloring
rule create_graph_coloring:
    input:
        ref_graph = "output/eval/par1_var/graphs/T2T.{num_samples}samples.{mapq}.{chrom}.L{min_var_size}kbp.gfa",
        seq_classes = "references_derived/T2T.chrY-seq-classes.tsv",
        dump = "references_derived/sample_annotations.json"
    output:
        table = "output/eval/par1_var/graphs/T2T.{num_samples}samples.{mapq}.{chrom}.L{min_var_size}kbp.annotations.csv"
    run:
        import pandas as pd
        import json

        seq_classes = pd.read_csv(
            input.seq_classes, sep="\t", header=0,
            usecols=["chrom", "start", "end", "name"]
        )

        with open(input.dump, "r") as dump:
            sample_infos = json.load(dump)

        node_infos = []
        with open(input.ref_graph, "r") as gfa:
            for line in gfa:
                if not line.startswith("S"):
                    break
                columns = line.strip().split()
                node_id = columns[1]
                sample_id = columns[4]
                assert sample_id.startswith("SN")
                sample_id = sample_id.split(":")[-1]
                if sample_id == "chrY":
                    sample_id = "T2T.Y"
                elif sample_id == "chrX":
                    sample_id = "T2T.X"
                else:
                    contig_infos = sample_id.split(".")
                    if contig_infos[0] == "chrY":
                        sample_id = contig_infos[-1] + ".Y"
                    elif contig_infos[0] == "chrX":
                        sample_id = contig_infos[-1] + ".X"
                    else:
                        raise ValueError(contig_infos)
                node_length = len(columns[2])
                if sample_id == "T2T.Y":
                    # that's the start coordinate
                    coord_offset = int(columns[5].split(":")[-1])
                    begin_idx = seq_classes.index[(seq_classes["start"] <= coord_offset)].max()
                    end_idx = seq_classes.index[(coord_offset + node_length <= seq_classes["end"])].min()
                    begin = seq_classes.at[begin_idx, "name"]
                    end = seq_classes.at[end_idx, "name"]
                    if begin == end:
                        seqclass = begin
                    else:
                        seqclass = f"{begin}-->>{end}"
                    start = coord_offset
                    end = coord_offset + node_length
                else:
                    seqclass = "non-ref"
                    start = -1
                    end = -1
                this_sample = dict(sample_infos[sample_id])
                this_sample["node"] = node_id
                this_sample["seqclass"] = seqclass
                this_sample["ref_start"] = start
                this_sample["ref_end"] = end
                node_infos.append(this_sample)
        df = pd.DataFrame.from_records(
            node_infos, columns=[
                "node", "color", "sample", "haplogroup",
                "population", "super_pop", "seqclass",
                "ref_start", "ref_end"
            ]
        )
        df.sort_values("node", inplace=True)
        df.to_csv(output.table, header=True, index=False)
    # END OF RUN BLOCK


rule run_all_graph_builds:
    input:
        gfa = expand(
            "output/eval/par1_var/graphs/T2T.{num_samples}samples.na.{chrom}.L{min_var_size}kbp.annotations.csv",
            num_samples=[3, 2, 1],
            chrom=["chrY", "chrXY"],
            min_var_size=[1, 5, 10]
        ),
