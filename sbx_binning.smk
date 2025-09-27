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
    threads: 8
    shell:
        """
        if [ -s {input.contigs} ]; then
            mkdir -p {output}
            metabat2 -i {input.contigs} -a {input.depth} -o {output}/{wildcards.sample} &> {log}
        else
            touch {output}/{wildcards.sample}.dummy
        fi
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
# Refinement (MAGScoT)
# ----------------------------
rule refine_bins:
    input:
        metabat2="bins/{sample}/metabat2",
        vamb="bins/{sample}/vamb"
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
        if [ -d {input.metabat2} ] || [ -d {input.vamb} ]; then
            mkdir -p $(dirname {output.bins})
            MAGScoT \
              --metabat {input.metabat2} \
              --vamb {input.vamb} \
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
