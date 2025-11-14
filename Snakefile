# Snakefile — cfNIPT pipelinREF_FASTAe
# Bowtie2 → Picard → SeqFF → idxstats → Wisecondor → WisecondorX → Report
# Usage:
#   snakemake --use-conda --conda-frontend conda --cores 15 --rerun-incomplete --keep-going
# or with your Slurm profile:
#   snakemake --profile resources/slurm --rerun-incomplete --keep-going

import os, glob, re

# ====================== CONFIG (edit if needed) ===============================
RUN_ID   = os.environ.get("RUN_ID", "*")   # include runs that have rawdata/<RUN>/.complete
THREADS  = 15
ROOT     = workflow.basedir

BWT2_PREFIX = f"{ROOT}/resources/genomeReferences/b37/bowtie2/human_g1k_v37_decoy"
REF_FASTA   = f"{ROOT}/resources/genomeReferences/b37/bowtie2/human_g1k_v37_decoy.fasta"
SEQFF_DIR   = f"{ROOT}/resources/software/easy-to-use-seqff"
PICARD_SIF  = f"{ROOT}/resources/software/picard-3.4.0--hdfd78af_0.sif"
WXC_SIF     = f"{ROOT}/resources/software/wisecondorx_1.2.9--pyhdfd78af_0.sif"
WXC_REF     = f"{ROOT}/resources/wisecondorReferences/wisecondorX_reference.npz"
WC_SIF      = f"{ROOT}/resources/software/nipt-wisecondor5.sif"
WC_REF      = f"{ROOT}/resources/wisecondorReferences/wisecondor_reference.npz"
REPORT_RMD  = f"{ROOT}/resources/reportgeneration.Rmd"
LOG_DIR     = f"{ROOT}/resources/logs"

# Conda envs
ENV_BWA    = "resources/envs/bwa.yml"            # bowtie2 + samtools
ENV_SEQFF  = "resources/envs/seqff-env.yaml"     # R + seqff deps
ENV_WC     = "resources/envs/bwa.yml"            # python2 + wisecondor (kept as-is)
ENV_R      = "resources/envs/r_env.yaml"         # rmarkdown + pkgs
ENV_FASTQC = "resources/envs/fastqc.yaml"         # fastqc

# Prevent locale spam from R
shell.prefix("set -euo pipefail; export LANG=C.UTF-8 LC_ALL=C.UTF-8; ")

# =================== DEFAULT TARGET FIRST (no wildcards) ======================
FINAL = []   # will be populated after sample discovery

rule all:
    input:
        lambda wc: FINAL
# ==============================================================================




# ==================== SIMPLE MULTI-SAMPLE + SINGLE-FILE DISCOVERY ============
import os, glob, re

fastq_files = sorted(glob.glob("rawdata/**/*.fastq.gz", recursive=True))
samples_ready = {}

for f in fastq_files:
    folder = os.path.dirname(f)
    base = os.path.basename(f)

    # --- CASE 1: Multi-lane (L001–L004) ---
    if base.endswith("_L001_R1_001.fastq.gz"):
        prefix = f[:-len("_L001_R1_001.fastq.gz")]
        lanes = [f"{prefix}_L00{i}_R1_001.fastq.gz" for i in range(1, 5)]
        if all(os.path.exists(p) for p in lanes):
            samples_ready[prefix] = folder

    # --- CASE 2: Any single FASTQ (no L00x pattern) ---
    elif not re.search(r"_L00[1-9]_R1_001\.fastq\.gz$", base):
        # Accept anything ending in .fastq.gz, even *_copy, *_rep, etc.
        prefix = re.sub(r"(?:_R1(?:_001)?)?\.fastq\.gz$", "", f)
        samples_ready[prefix] = folder

if not samples_ready:
    print("[INFO] No complete or single FASTQs found under rawdata/")
else:
    print("[INFO] Found samples ready for processing:")
    for s, folder in samples_ready.items():
        print(f"   - {os.path.basename(s)}  ({folder})")

# Collect all FASTQs (either 4 lanes or 1 file per sample)
_all_fastqs = []
for s in samples_ready:
    lanes = [f"{s}_L00{i}_R1_001.fastq.gz" for i in range(1, 5)]
    if all(os.path.exists(p) for p in lanes):
        _all_fastqs.extend(lanes)
    else:
        # Add first existing single-file variant
        for pat in (f"{s}_R1_001.fastq.gz", f"{s}_R1.fastq.gz", f"{s}.fastq.gz"):
            if os.path.exists(pat):
                _all_fastqs.append(pat)
                break

# Derive unique sample names
SAMPLES = sorted({
    re.sub(r'(?:_L00[1-9])?(?:_R1(?:_001)?)?$', '', os.path.basename(p)[:-9])
    for p in _all_fastqs
})

# Compatibility for concat_lanes logic
RUN_DIRS = sorted(set(samples_ready.values()))
# ==============================================================================








# ======================= SIMPLE L001–L004 CONCAT ==============================
def _simple_lane_paths(sample):
    """
    If any {sample}_L001_R1_001.fastq.gz exists in a run, concatenate
    any of L001..L004 present in THAT SAME run (order preserved).
    Else fall back to a single-file input in priority: _R1_001, _R1, bare .fastq.gz.
    """
    # prefer a run that has at least one lane file
    for run in RUN_DIRS:
        lanes = [f"{run}/{sample}_L00{i}_R1_001.fastq.gz" for i in range(1, 5)]
        present = [p for p in lanes if os.path.exists(p)]
        if present:
            # ensure order L001..L004 by the lanes list order
            return present
    # fallback: single-file patterns, first match wins across runs
    for run in RUN_DIRS:
        for pat in (f"{sample}_R1_001.fastq.gz",
                    f"{sample}_R1.fastq.gz",
                    f"{sample}.fastq.gz"):
            p = f"{run}/{pat}"
            if os.path.exists(p):
                return [p]
    searched = ", ".join(RUN_DIRS) if RUN_DIRS else "(no eligible runs)"
    raise ValueError(f"No FASTQs found for sample '{sample}' in: {searched}")

def lanes_or_single(wildcards):
    return _simple_lane_paths(wildcards.sample)

rule concat_lanes:
    input:
        lanes_or_single
    output:
        "results/{sample}/raw/{sample}.fastq.gz"
    threads: 1
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        cat {input} > "{output}"
        """
# ==============================================================================


# ====================== BUILD FINAL TARGETS (per sample) ======================
for s in SAMPLES:
    FINAL += [
        f"results/{s}/raw/{s}.fastq.gz",                          # concatenated/copy FASTQ
        f"results/{s}/alignments/{s}.bowtie2-b37.bam",
        f"results/{s}/alignments/{s}.bowtie2-b37.bam.bai",
        f"results/{s}/intermediateOutput/{s}.gc_bias_metrics.txt",
        f"results/{s}/intermediateOutput/{s}.gc_bias_metrics.pdf",
        f"results/{s}/intermediateOutput/{s}.gc_summary_metrics.txt",
        f"results/{s}/intermediateOutput/{s}.seqff.txt",
        f"results/{s}/intermediateOutput/{s}.bowtie2-b37-idxstats.txt",
        f"results/{s}/alignments/{s}.wisecondor-bowtie2_b37.npz",
        f"results/{s}/alignments/{s}.wisecondor-bowtie2_b37-out.npz",
        f"results/{s}/intermediateOutput/{s}.wisecondor-bowtie2_b37-report.txt",
        f"results/{s}/intermediateOutput/{s}-wisecondor-bowtie2_b37_z.pdf",
        f"results/{s}/intermediateOutput/{s}.wisecondor-bowtie2_b37_cwz.csv",
        f"results/{s}/intermediateOutput/{s}_wisecondorx.npz",
        f"results/{s}/intermediateOutput/{s}.wisecondorx-b37_gender.txt",
        f"results/{s}/intermediateOutput/{s}_report.pdf",
        f"PDF/{s}_report.pdf",
        f"results/{s}/raw/{s}_fastqc.html",
    ]
# ==============================================================================


# ========================= ALIGN → SORTED BAM =================================
rule bowtie2_align_sort:
    input:
        fastq = "results/{sample}/raw/{sample}.fastq.gz"   # use concatenated/copy
    output:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    params:
        idx = BWT2_PREFIX
    threads: THREADS
    resources:
        mem_mb = 36000
    conda: ENV_BWA
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        bowtie2 --local --fast-local -x "{params.idx}" --threads {threads} -U "{input.fastq}" \
          | samtools sort -@ {threads} -m 2G -O BAM -o "{output}"
        """
# ==============================================================================


# =============================== BAM INDEX ====================================
rule samtools_index:
    input:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    output:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam.bai"
    threads: 4
    conda: ENV_BWA
    shell:
        r"""samtools index -@ {threads} "{input}" "{output}" """
# ==============================================================================


# ============================== PICARD GC BIAS ================================
rule picard_gc_bias:
    input:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    output:
        "results/{sample}/intermediateOutput/{sample}.gc_bias_metrics.txt",
        "results/{sample}/intermediateOutput/{sample}.gc_bias_metrics.pdf",
        "results/{sample}/intermediateOutput/{sample}.gc_summary_metrics.txt"
    params:
        fasta = REF_FASTA,
        sif   = PICARD_SIF
    shell:
        r"""
        mkdir -p "$(dirname {output[0]})"
        apptainer exec "{params.sif}" \
          picard CollectGcBiasMetrics \
            -I "{input}" \
            -O "{output[0]}" \
            -CHART "{output[1]}" \
            -S "{output[2]}" \
            -R "{params.fasta}"
        """
# ==============================================================================


# ================================ SeqFF =======================================
rule seqff:
    input:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    output:
        "results/{sample}/intermediateOutput/{sample}.seqff.txt"
    params:
        tool = SEQFF_DIR,
        root = ROOT
    conda: ENV_SEQFF
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        cd "{params.tool}"
        Rscript seqff.r -f "{params.root}/{input}" -o "{params.root}/{output}"
        """
# ==============================================================================


# ============================== samtools idxstats =============================
rule samtools_idxstats:
    input:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    output:
        "results/{sample}/intermediateOutput/{sample}.bowtie2-b37-idxstats.txt"
    conda: ENV_BWA
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        samtools idxstats "{input}" > "{output}"
        """
# ==============================================================================


# ============================== WISECONDOR (py2) ==============================
rule wisecondor_convert_b37:
    input:
        "results/{sample}/alignments/{sample}.bowtie2-b37.bam"
    output:
        "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37.npz"
    params:
        sif = WC_SIF,
        binsize = 50000
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer -q run --cleanenv "{params.sif}" \
          convert "{input}" "{output}" -binsize {params.binsize}
        """

rule wisecondor_test_b37:
    input:
        "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37.npz"
    output:
        "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37-out.npz"
    params:
        sif = WC_SIF,
        ref = WC_REF,
        minz = 9.9
    shell:
        r"""
        apptainer -q run --cleanenv "{params.sif}" \
          test "{input}"  "{output}" "{params.ref}" -minzscore {params.minz}
        """

rule wisecondor_report_b37:
    input:
        npz    = "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37.npz",
        npzout = "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37-out.npz"
    output:
        "results/{sample}/intermediateOutput/{sample}.wisecondor-bowtie2_b37-report.txt"
    params:
        sif = WC_SIF
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer -q run --cleanenv "{params.sif}" \
          report "{input.npz}" "{input.npzout}" > "{output}"
        """

rule wisecondor_plot_b37:
    input:
        "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37-out.npz"
    output:
        "results/{sample}/intermediateOutput/{sample}-wisecondor-bowtie2_b37_z.pdf"
    params:
        sif = WC_SIF,
        prefix = "results/{sample}/intermediateOutput/{sample}-wisecondor-bowtie2_b37"
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer -q run --cleanenv "{params.sif}" \
          plot "{input}" "{params.prefix}"
        """

rule wisecondor_extract_cwz_b37:
    input:
        "results/{sample}/alignments/{sample}.wisecondor-bowtie2_b37-out.npz"
    output:
        "results/{sample}/intermediateOutput/{sample}.wisecondor-bowtie2_b37_cwz.csv"
    params:
        sif = WC_SIF
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer -q exec --cleanenv "{params.sif}" \
          wisecondor-cwz "{input}" "{output}"
        """
# ==============================================================================


# =============================== WISECONDORX ==================================
rule wisecondorx_convert:
    input:
        bam = "results/{sample}/alignments/{sample}.bowtie2-b37.bam",
        bai = "results/{sample}/alignments/{sample}.bowtie2-b37.bam.bai"
    output:
        "results/{sample}/intermediateOutput/{sample}_wisecondorx.npz"
    params:
        sif    = WXC_SIF,
        prefix = "results/{sample}/intermediateOutput/{sample}_wisecondorx"
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer exec "{params.sif}" \
          WisecondorX convert "{input.bam}" "{params.prefix}"
        """

rule wisecondorx_predict:
    input:
        "results/{sample}/intermediateOutput/{sample}_wisecondorx.npz"
    output:
        bins        = "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_bins.bed",
        segments    = "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_segments.bed",
        aberrations = "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_aberrations.bed",
        statistics  = "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_statistics.txt",
        plots       = directory("results/{sample}/intermediateOutput/{sample}.wisecondorx-b37.plots")
    params:
        sif = WXC_SIF,
        ref = WXC_REF,
        prefix = "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37"
    shell:
        r"""
        mkdir -p "$(dirname {output.bins})"
        apptainer exec "{params.sif}" \
          WisecondorX predict "{input}" "{params.ref}" "{params.prefix}" --bed --plot
        """

rule wisecondorx_gender:
    input:
        "results/{sample}/intermediateOutput/{sample}_wisecondorx.npz"
    output:
        "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_gender.txt"
    params:
        sif = WXC_SIF,
        ref = WXC_REF
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        apptainer exec "{params.sif}" \
          WisecondorX gender "{input}" "{params.ref}" > "{output}"
        """
# ==============================================================================


# =============================== REPORT (Rmd) =================================
rule report_pdf:
    input:
        "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_gender.txt",
        "results/{sample}/intermediateOutput/{sample}.wisecondorx-b37_statistics.txt"
    output:
        "results/{sample}/intermediateOutput/{sample}_report.pdf"
    params:
        rmd = REPORT_RMD
    conda: ENV_R
    shell:
        r"""
        mkdir -p "$(dirname {output})" "results/{wildcards.sample}/rmd_tmp" "{LOG_DIR}"
        Rscript -e "rmarkdown::render('{params.rmd}', \
          output_format='pdf_document', \
          output_file=basename('{output}'), \
          output_dir=dirname('{output}'), \
          params=list(fileid='{wildcards.sample}'), \
          intermediates_dir='results/{wildcards.sample}/rmd_tmp')" \
          "{wildcards.sample}"
        """
# ==============================================================================


# =============================== COPY REPORT ================================== 
rule copy_pdf:
    input:
        "results/{sample}/intermediateOutput/{sample}_report.pdf"
    output:
        "PDF/{sample}_report.pdf"
    threads: 1
    resources:
        mem_mb=512
    shell:
        r"""
        cp "{input}" "{output}"
        """
# ==============================================================================


# ================================== FASTQC ====================================
rule fastqc_stats:
    input:
        "results/{sample}/raw/{sample}.fastq.gz"
    output:
        "results/{sample}/raw/{sample}_fastqc.html"
    conda: ENV_FASTQC
    threads: 2
    resources:
        mem_mb=1024
    shell:
        r"""
        fastqc "{input}"; echo "{output}"
        """
# ==============================================================================


