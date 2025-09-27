try:
    SBX_ASSEMBLY_VERSION = get_ext_version("sbx_assembly")
except NameError:
    SBX_ASSEMBLY_VERSION = "0.0.0"

# Collect all outputs in one rule
rule all_binning:
    input:
        expand("bins/{sample}/refined/{sample}.refined_bins.fa", sample=Samples),
        expand("qc/mags/{sample}.checkm2.tsv", sample=Samples),

# ----------------------------
# Depth summarization
# ----------------------------
rule contig_depth_summary:
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        bams=expand(ASSEMBLY_FP / "coverage" / "samtools" / "{other}.sorted.bam",
                    other=Samples.keys())
    output:
        depth=ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv",
    benchmark:
        BENCHMARK_FP / "contig_depth_summary_{sample}.tsv"
    log:
        LOG_FP / "contig_depth_summary_{sample}.log",
    conda:
        "envs/sbx_binning.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
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
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        depth=ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv",
    output:
        directory("bins/{sample}/metabat2"),
    benchmark:
        BENCHMARK_FP / "binning_metabat2_{sample}.tsv"
    log:
        LOG_FP / "binning_metabat2_{sample}.log",
    conda:
        "envs/sbx_binning.yml"
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
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        depth=ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv",
    output:
        directory("bins/{sample}/vamb"),
    benchmark:
        BENCHMARK_FP / "binning_vamb_{sample}.tsv"
    log:
        LOG_FP / "binning_vamb_{sample}.log",
    conda:
        "envs/sbx_binning.yml"
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


# ----------------------------
# Make scaffolds2bin.tsv from MetaBAT2 bins
# ----------------------------
rule metabat2_scaffolds2bin:
    input:
        bins_dir="bins/{sample}/metabat2"
    output:
        tsv="bins/{sample}/metabat2_scaffolds2bin.tsv"
    benchmark:
        BENCHMARK_FP / "metabat2_scaffolds2bin_{sample}.tsv"
    log:
        LOG_FP / "metabat2_scaffolds2bin_{sample}.log",
    conda:
        "envs/sbx_binning.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    shell:
        """
        if compgen -G "{input.bins_dir}/*.fa" > /dev/null; then
            : > {output.tsv}
            for F in {input.bins_dir}/*.fa; do
              BIN=$(basename $F .fa)
              awk -v b=$BIN '/^>/{gsub(/^>/, "", $0); split($0,a,/[\t ]/); print a[1]"\t"b}' $F >> {output.tsv}
            done &> {log}
        else
            touch {output.tsv}
        fi
        """

# ----------------------------
# Convert VAMB clusters.tsv to scaffolds2bin.tsv
# ----------------------------
rule vamb_scaffolds2bin:
    input:
        clusters="bins/{sample}/vamb/clusters.tsv"
    output:
        tsv="bins/{sample}/vamb_scaffolds2bin.tsv"
    benchmark:
        BENCHMARK_FP / "vamb_scaffolds2bin_{sample}.tsv"
    log:
        LOG_FP / "vamb_scaffolds2bin_{sample}.log",
    conda:
        "envs/sbx_binning.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    script:
        "scripts/vamb_clusters_to_bins.py"


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
        "envs/sbx_binning.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-binning"
    threads: 8
    shell:
        """
        if [ -s {input.metabat2_map} ] || [ -s {input.vamb_map} ]; then
            mkdir -p $(dirname {output.bins})
            MAGScoT \
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
        "envs/sbx_binning.yml"
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
