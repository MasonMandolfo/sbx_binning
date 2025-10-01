try:
    SBX_ASSEMBLY_VERSION = get_ext_version("sbx_assembly")
except NameError:
    SBX_ASSEMBLY_VERSION = "0.0.0"

# ----------------------------
# Top-level target for binning workflow
# ----------------------------
# ----------------------------
# Top-level target for binning workflow
# ----------------------------
rule all_binning:
    input:
        expand("bins/{sample}/refined/{sample}.refined_bins.fa", sample=Samples),
        expand("qc/mags/{sample}.checkm2.tsv", sample=Samples),


import csv

def load_groups(csv_fp):
    groups = {}
    with open(csv_fp) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            group_id, *samples = row
            groups[group_id] = samples
    return groups

group_csv = Cfg["sbx_binning"].get("group_csv")
if not group_csv:
    raise ValueError("Please specify sbx_binning.group_csv in sunbeam_config.yml")

Groups = load_groups(group_csv)

# Build sample→group index
Sample2Group = {s: g for g, samples in Groups.items() for s in samples}

def group_members(sample):
    """Return the list of samples in the same group as `sample` (anchor)."""
    g = Sample2Group.get(sample)
    if g is None:
        raise ValueError(f"Sample {sample} not found in group_csv {group_csv}")
    return Groups[g]

rule crossmap_sort:
    input:
        contig = ASSEMBLY_FP / "contigs" / "{anchor}-contigs.fa",
        reads  = expand(QC_FP / "decontam" / "{{other}}_{rp}.fastq.gz", rp=Pairs)
    output:
        bam = ASSEMBLY_FP / "coverage" / "samtools_cross" / "{anchor}__{other}.sorted.bam"
    log:
        LOG_FP / "crossmap_sort_{anchor}__{other}.log"
    threads: Cfg["sbx_assembly"]["threads"]
    conda:
        "envs/sbx_binning_env.yml"   
    shell:
        r"""
        # align then sort
        mkdir -p $(dirname {output.bam})
        minimap2 -ax sr -t {threads} {input.contig} {input.reads} \
          2> {log} \
        | samtools sort -@ {threads} -o {output.bam} -  >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """

def cross_bams_for_anchor(wc):
    """All group members' BAMs mapped to the anchor's contigs."""
    return [
        ASSEMBLY_FP / "coverage" / "samtools_cross" / f"{wc.anchor}__{other}.sorted.bam"
        for other in group_members(wc.anchor)
    ]

rule summarize_bam_depths:
    input:
        contigs = ASSEMBLY_FP / "contigs" / "{anchor}-contigs.fa",
        bams    = cross_bams_for_anchor,
    output:
        depth = ASSEMBLY_FP / "coverage" / "depth" / "{anchor}.contig_depth.tsv",
    log:
        LOG_FP / "summarize_bam_depths_{anchor}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 4
    shell:
        """
        if [ -s {input.contigs} ]; then
            mkdir -p $(dirname {output.depth})
            jgi_summarize_bam_contig_depths \
                --outputDepth {output.depth} \
                {input.bams} &> {log}
        else
            touch {output.depth}
        fi
        """


# ----------------------------
# MetaBAT2
# ----------------------------
rule binning_metabat2:
    input:
        contigs =   ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        depth   =   ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv"
    output:
        directory("bins/{sample}/metabat2"),
    benchmark:
        BENCHMARK_FP / "binning_metabat2_{sample}.tsv"
    log:
        LOG_FP / "binning_metabat2_{sample}.log",
    conda:
        "envs/sbx_binning_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    params:
        min_contig=Cfg["sbx_binning"]["min_contig"]
    shell:
        """
        metabat2 -i {input.contigs} -a {input.depth} \
                -m {params.min_contig} \
                -o {output}/{wildcards.sample} &> {log}
        """


# ----------------------------
# VAMB
# ----------------------------
rule binning_vamb:
    input:
        contigs =   ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        depth   =   ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv"
    output:
        clusters="bins/{sample}/vamb/clusters.tsv"
    benchmark:
        BENCHMARK_FP / "binning_vamb_{sample}.tsv"
    log:
        LOG_FP / "binning_vamb_{sample}.log",
    conda:
        "envs/sbx_binning_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    threads: 16
    shell:
        """
        if [ -s {input.contigs} ]; then
            mkdir -p {output}
            vamb --outdir {output} \
                 --fasta {input.contigs} \
                 --jgi {input.depth} \
                 --minfasta 200000 &> {log}
        else
            touch {output}/vamb.dummy
        fi
        """

rule scaffolds2bin:
    input:
        # Use lambda to select correct input based on tool
        lambda wildcards: (
            f"bins/{wildcards.sample}/metabat2"
            if wildcards.tool == "metabat2"
            else f"bins/{wildcards.sample}/vamb/clusters.tsv"
        )
    output:
        tsv="bins/{sample}/{tool}_scaffolds2bin.tsv"
    log:
        LOG_FP / "scaffolds2bin_{sample}_{tool}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 1
    run:
        import os, subprocess
        if wildcards.tool == "metabat2":
            bins_dir = input[0]
            with open(output.tsv, "w") as out, open(log, "w") as lg:
                if any(fn.endswith(".fa") for fn in os.listdir(bins_dir)):
                    for fn in os.listdir(bins_dir):
                        if fn.endswith(".fa"):
                            path = os.path.join(bins_dir, fn)
                            BIN = os.path.splitext(fn)[0]
                            cmd = f"awk -v b={BIN} '/^>/' {path}"
                            res = subprocess.check_output(cmd, shell=True, text=True)
                            for line in res.strip().splitlines():
                                contig = line.lstrip(">").split()[0]
                                out.write(f"{contig}\t{BIN}\n")
        elif wildcards.tool == "vamb":
            shell(f"python scripts/vamb_clusters_to_bins.py {input} {output.tsv} > {log} 2>&1")


# ----------------------------
# Refinement (MAGScoT)
# ----------------------------
rule refine_bins:
    input:
        metabat2_bins="bins/{sample}/metabat2",
        metabat2_map="bins/{sample}/metabat2_scaffolds2bin.tsv",
        vamb_map="bins/{sample}/vamb_scaffolds2bin.tsv"
    output:
        bins="bins/{sample}/refined/{sample}.refined_bins.fa",
        summary="bins/{sample}/refined/{sample}.refined_summary.tsv"
    benchmark:
        BENCHMARK_FP / "refine_bins_{sample}.tsv"
    log:
        LOG_FP / "refine_bins_{sample}.log",
    conda:
        "envs/sbx_binning_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    threads: 8
    shell:
        """
        if [ -s {input.metabat2_map} ] || [ -s {input.vamb_map} ]; then
            mkdir -p $(dirname {output.bins})
             {workflow.basedir}/scripts/run_magscot.sh \
              --metabat {input.metabat2_map} \
              --vamb {input.vamb_map} \
              --outdir $(dirname {output.bins}) \
              --prefix {wildcards.sample} &> {log}
        else
            touch {output.bins} {output.summary}
        fi
        """


# ----------------------------
# QC (CheckM2)
# ----------------------------
rule qc_bins:
    input:
        bins="bins/{sample}/refined/{sample}.refined_bins.fa"
    output:
        tsv="qc/mags/{sample}.checkm2.tsv"
    benchmark:
        BENCHMARK_FP / "qc_bins_{sample}.tsv"
    log:
        LOG_FP / "qc_bins_{sample}.log",
    conda:
        "envs/sbx_checkm2_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    threads: 8
    shell:
        """
        if [ -s {input.bins} ]; then
            mkdir -p qc/mags
            checkm2 predict \
              --input {input.bins} \
              --output qc/mags/{wildcards.sample}.checkm2 \
              --threads {threads} &> {log}
            cp qc/mags/{wildcards.sample}.checkm2/quality_report.tsv {output.tsv}
        else
            touch {output.tsv}
        fi
        """
