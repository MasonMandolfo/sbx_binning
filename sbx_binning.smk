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

rule depth_to_metabat2:
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

rule depth_to_vamb:
    input:
        depth = ASSEMBLY_FP / "coverage" / "depth" / "{sample}.contig_depth.tsv"
    output:
        vambdepth = ASSEMBLY_FP / "coverage" / "depth" / "{sample}.vamb_depth.tsv"
    log:
        LOG_FP / "depth_to_vamb_{sample}.log"
    run:
        with open(input.depth) as fin, open(output.vambdepth, "w") as fout:
            header = fin.readline().strip().split("\t")
            # Keep only the contigName and the BAM mean-depth columns (skip contigLen, totalAvgDepth, and -var)
            keep_idxs = [0] + [i for i, h in enumerate(header) if h.endswith(".sorted.bam") and not h.endswith("-var")]
            # Write new header: "contigname" + clean sample IDs
            samples = [h.split("__")[-1].replace(".sorted.bam", "") for h in header if h.endswith(".sorted.bam")]
            fout.write("contigname\t" + "\t".join(samples) + "\n")
            for line in fin:
                parts = line.strip().split("\t")
                fout.write("\t".join(parts[i] for i in keep_idxs) + "\n")



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
        depth   =   ASSEMBLY_FP / "coverage" / "depth" / "{sample}.vamb_depth.tsv"
    output:
        vambdir = directory("bins/{sample}/vamb")
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
            vamb bin default\
                 --outdir {output.vambdir}  \
                 --fasta {input.contigs} \
                 --abundance_tsv {input.depth} \
                 --minfasta 200000 &> {log}
        else
            touch {output.vambdir}/vamb.dummy
        fi
        """

rule scaffolds2bin:
    input:
        lambda wildcards: (
            f"bins/{wildcards.sample}/metabat2"
            if wildcards.tool == "metabat2"
            else directory(f"bins/{wildcards.sample}/vamb")
        )
    output:
        tsv = "bins/{sample}/{tool}_scaffolds2bin.tsv"
    log:
        LOG_FP / "scaffolds2bin_{sample}_{tool}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 1
    shell:
        r"""
        if [ "{wildcards.tool}" = "metabat2" ]; then
            bins_dir={input}
            if ls ${{bins_dir}}/*.fa 1> /dev/null 2>&1; then
                for f in ${{bins_dir}}/*.fa; do
                    bin=$(basename $f .fa)
                    awk -v b=$bin '/^>/ {{gsub(/^>/, "", $1); print b "\t" $1 "\tmetabat2"}}' $f
                done > {output.tsv} 2> {log}
            else
                touch {output.tsv}
            fi
        elif [ "{wildcards.tool}" = "vamb" ]; then
            awk 'NR>1 {{print $2 "\t" $1 "\tvamb"}}' {input} > {output.tsv} 2> {log}
        else
            echo "Unknown binning tool: {wildcards.tool}" >&2
            exit 1
        fi
        """

rule concat_bin_mappings:
    input:
        met = "bins/{sample}/metabat2_scaffolds2bin.tsv",
        vam = "bins/{sample}/vamb_scaffolds2bin.tsv"
    output:
        concat = "bins/{sample}/contigs_to_bin.tsv"
    log:
        LOG_FP / "concat_bin_mappings_{sample}.log"
    shell:
        r"""
        set -euo pipefail
        # Ensure proper 3-column format: BinID ContigID Tool
        cat {input.met} {input.vam} > {output.concat} 2> {log}
        """

rule prodigal_orf:
    input:
        contigs = ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa"
    output:
        faa="qc/{sample}/prodigal/{sample}.faa",
        ffn="qc/{sample}/prodigal/{sample}.ffn"
    log:
        LOG_FP / "prodigal_{sample}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 8
    shell:
       r"""
        set -euo pipefail
        mkdir -p $(dirname {output.faa})
        tmpdir="$(mktemp -d)"

        # Split FASTA by records and run prodigal in parallel; concatenate chunk outputs
        cat {input.contigs} \
          | parallel -j {threads} --block 999k --recstart '>' --pipe \
              prodigal -p meta \
                       -a "$tmpdir"/{wildcards.sample}.{{#}}.faa \
                       -d "$tmpdir"/{wildcards.sample}.{{#}}.ffn \
                       -o /dev/null 2>> {log}

        cat "$tmpdir"/{wildcards.sample}.*.faa > {output.faa}
        cat "$tmpdir"/{wildcards.sample}.*.ffn > {output.ffn}
        rm -rf "$tmpdir"
        """

MAGScoT_fp = Cfg["magscot"]["MAGScoT_fp"]

rule hmmsearch_tigr:
    input:
        faa="qc/{sample}/prodigal/{sample}.faa"
    output:
        tbl = "qc/{sample}/hmm/{sample}.hmm.tigr.hit.out",
        out = "qc/{sample}/hmm/{sample}.hmm.tigr.out"
    log:
        LOG_FP / "hmmsearch_tigr_{sample}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl})
        hmmsearch -o {output.out} --tblout {output.tbl} \
          --noali --notextw --cut_nc --cpu {threads} \
          "{MAGScoT_fp}/hmm/gtdbtk_rel207_tigrfam.hmm" \
          {input.faa} &> {log}
        """

rule hmmsearch_pfam:
    input:
        faa="qc/{sample}/prodigal/{sample}.faa"
    output:
        tbl = "qc/{sample}/hmm/{sample}.hmm.pfam.hit.out",
        out = "qc/{sample}/hmm/{sample}.hmm.pfam.out"
    log:
        LOG_FP / "hmmsearch_pfam_{sample}.log"
    conda:
        "envs/sbx_binning_env.yml"
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl})
        hmmsearch -o {output.out} --tblout {output.tbl} \
          --noali --notextw --cut_nc --cpu {threads} \
          "{MAGScoT_fp}/hmm/gtdbtk_rel207_Pfam-A.hmm" \
          {input.faa} &> {log}
        """

rule combine_hmm_hits:
    input:
        tigr="qc/{sample}/hmm/{sample}.hmm.tigr.hit.out",
        pfam="qc/{sample}/hmm/{sample}.hmm.pfam.hit.out"
    output:
        hmm="qc/{sample}/hmm/{sample}.hmm"
    shell:
        r"""
        (grep -v "^#" {input.tigr} | awk '{{print $1"\t"$3"\t"$5}}'
        grep -v "^#" {input.pfam} | awk '{{print $1"\t"$4"\t"$5}}') > {output.hmm}
        """


# ----------------------------
# Refinement (MAGScoT)
# ----------------------------


rule run_magscot:
    input:
        contig_map = "bins/{sample}/contigs_to_bin.tsv",
        hmm = "qc/{sample}/hmm/{sample}.hmm"
    output:
        outdir = directory("qc/{sample}/refined/{sample}.magscot"),
        sentinel="qc/{sample}/refined/{sample}.magscot/done"
    log:
        LOG_FP / "magscot_{sample}.log"
    params:
        profile = Cfg.get("magscot", {}).get("profile", "bac120+ar53"),
        score_a = lambda wc: f"--a {Cfg['magscot']['score_a']}" if Cfg['magscot'].get('score_a') else "",
        score_b = lambda wc: f"--b {Cfg['magscot']['score_b']}" if Cfg['magscot'].get('score_b') else "",
        score_c = lambda wc: f"--c {Cfg['magscot']['score_c']}" if Cfg['magscot'].get('score_c') else "",
        max_cont = lambda wc: f"--max_cont {Cfg['magscot']['max_cont']}" if Cfg['magscot'].get('max_cont') else "",
        threshold = lambda wc: f"-t {Cfg['magscot']['threshold']}" if Cfg['magscot'].get('threshold') else "",
        skip_merge = lambda wc: "--skip_merge_bins" if Cfg['magscot'].get('skip_merge_bins', False) else "",
        magscot_fp = MAGScoT_fp
    conda:
        "envs/sbx_binning_env.yml"
    threads: Cfg["magscot"].get("threads", 4)
    shell:
        r"""
        set -euo pipefail
        Rscript {params.magscot_fp}/MAGScoT.R \
            -i {input.contig_map} \
            --hmm {input.hmm} \
            -p {params.profile} \
            -o {output.outdir} \
            {params.score_a} \
            {params.score_b} \
            {params.score_c} \
            {params.max_cont} \
            {params.threshold} \
            {params.skip_merge} \
            > {log} 2>&1
        touch {output.sentinel}
        """



# ----------------------------
# QC (CheckM2)
# ----------------------------
rule qc_bins:
    input:
        sentinel="qc/{sample}/refined/{sample}.magscot/done"
    output:
        tsv="qc/mags/{sample}.checkm2.tsv"
    log:
        LOG_FP / "qc_bins_{sample}.log"
    conda:
        "envs/sbx_checkm2_env.yml"
    threads: 8
    shell:
        r"""
        magscot_dir=$(dirname {input.sentinel})
        if compgen -G "$magscot_dir/refined_bins/*.fa" > /dev/null; then
            mkdir -p qc/mags
            checkm2 predict \
              --input $magscot_dir/refined_bins \
              --output qc/mags/{wildcards.sample}.checkm2 \
              --threads {threads} &> {log}
            cp qc/mags/{wildcards.sample}.checkm2/quality_report.tsv {output.tsv}
        else
            touch {output.tsv}
        fi
        """
