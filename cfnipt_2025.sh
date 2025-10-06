#!/bin/bash
#SBATCH --account KGA
#SBATCH -c 15
#SBATCH --mem 32g
#SBATCH --time 48:00:00
#SBATCH --job-name=cfnipt_2025
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err


# For --report-only mode: sbatch -c 2 --mem=8g ./cfnipt_2025.sh VALIDATION_150925 1028-13_S5_R1 --report-only
# Or for f in results/*; do sbatch -c 2 --mem 8g ./cfnipt_2025.sh VALIDATION_150925 $(basename $f) --report-only; done
# Or for multiple jobs:  for f in rawData2/*/*fastq.gz; do a=$(basename $f); sbatch ./cfnipt_2025.sh VALIDATION_150925 ${a%.fastq.gz}; done
set -Eeuo pipefail

########################
# Input argument
REFERENCE=analysisScripts/genomeReferences/b37/bowtie2/human_g1k_v37_decoy

RUN_ID="${1:?Usage: $0 RUN_ID SAMPLE [OUTDIR] [FASTQ_OVERRIDE]}"
FILEID="${2:?Usage: $0 RUN_ID SAMPLE [OUTDIR] [FASTQ_OVERRIDE]}"

# Optional overrides (keep defaults)
RESULTSFOLDER="${4:-results}"
MODE="${3:-}" 
FASTQ_OVERRIDE="${5:-}"

# Threads: prefer SLURM value if available
CORES="${SLURM_CPUS_PER_TASK:-15}"

# FASTQ detection (minimal + robust)
if [[ -n "$FASTQ_OVERRIDE" ]]; then
  FASTQ_IN="$FASTQ_OVERRIDE"
else
  base="rawData/$RUN_ID"
  if   [[ -f "$base/${FILEID}.fastq.gz"     ]]; then FASTQ_IN="$base/${FILEID}.fastq.gz"
  elif [[ -f "$base/${FILEID}_R1.fastq.gz"  ]]; then FASTQ_IN="$base/${FILEID}_R1.fastq.gz"
  elif [[ -f "$base/${FILEID}_R1_001.fastq.gz" ]]; then FASTQ_IN="$base/${FILEID}_R1_001.fastq.gz"
  else
    echo "ERROR: FASTQ not found for $FILEID in $base" >&2
    exit 2
  fi
fi

BAMFILE="$RESULTSFOLDER/$FILEID/alignments/${FILEID}.bowtie2-b37.bam"
RES_INTERMEDIATE="$RESULTSFOLDER/$FILEID/intermediateOutput"
########################

source "/home/$USER/miniforge3/etc/profile.d/conda.sh"

if ! [[ "$MODE" == "--report-only" ]]; then

conda activate bwa

mkdir -p "$RES_INTERMEDIATE" "$RESULTSFOLDER/$FILEID/alignments"

#######################
# BOWTIE2 ALIGNMENT
#######################
bowtie2 --local --fast-local -x "$REFERENCE" --threads "$CORES" -U "$FASTQ_IN" \
  | samtools sort -@ "$CORES" -m 2G -O BAM - > "$BAMFILE"

#######################
# INDEX
#######################
samtools index "$BAMFILE"

#######################
# PICARD
#######################
# Using container from: apptainer pull picard-3.4.0--hdfd78af_0.sif docker://quay.io/biocontainers/picard:3.4.0--hdfd78af_0
# Installed in folder analysisScripts/software/
apptainer exec analysisScripts/software/picard-3.4.0--hdfd78af_0.sif \
  picard CollectGcBiasMetrics \
    -I "$BAMFILE" \
    -O "$RES_INTERMEDIATE/${FILEID}.gc_bias_metrics.txt" \
    -CHART "$RES_INTERMEDIATE/${FILEID}.gc_bias_metrics.pdf" \
    -S "$RES_INTERMEDIATE/${FILEID}.gc_summary_metrics.txt" \
    -R /faststorage/project/KGA/cfNIPT_Wisecondor/analysisScripts/genomeReferences/b37/bowtie2/human_g1k_v37_decoy.fasta

conda deactivate
conda activate seqff-env

#######################
# SEQFF
#######################
# Directory is changed because of dependencies in the R code
cd analysisScripts/software/easy-to-use-seqff/
Rscript seqff.r -f "../../../$BAMFILE" -o "../../../$RES_INTERMEDIATE/${FILEID}.seqff.txt"
cd ../../../

### Using container from: apptainer pull wisecondorx_1.2.9--pyhdfd78af_0.sif docker://quay.io/biocontainers/wisecondorx:1.2.9--pyhdfd78af_0

conda deactivate
conda activate bwa

#######################
# SAMTOOLS IDXSTAT
#######################
samtools idxstats "$BAMFILE" > "$RES_INTERMEDIATE/${FILEID}.bowtie2-b37-idxstats.txt"

#######################
# WISECONDOR
#######################
snakemake -s analysisScripts/scripts/snakemake-nipt-v2.1-wisecondor -p \
  --config runID="$RUN_ID" sample="$FILEID" outdir="$RESULTSFOLDER" \
  --use-conda --rerun-incomplete --cores "$CORES"

#######################
# WISECONDORX
#######################
apptainer exec analysisScripts/software/wisecondorx_1.2.9--pyhdfd78af_0.sif \
  WisecondorX convert "$BAMFILE" "$RES_INTERMEDIATE/${FILEID}_wisecondorx"

apptainer exec analysisScripts/software/wisecondorx_1.2.9--pyhdfd78af_0.sif \
  WisecondorX predict "$RES_INTERMEDIATE/${FILEID}_wisecondorx.npz" \
    analysisScripts/wisecondorReferences/wisecondorX_reference.npz \
    "$RES_INTERMEDIATE/${FILEID}.wisecondorx-b37" --bed --plot

apptainer exec analysisScripts/software/wisecondorx_1.2.9--pyhdfd78af_0.sif \
  WisecondorX gender "$RES_INTERMEDIATE/${FILEID}_wisecondorx.npz" \
    analysisScripts/wisecondorReferences/wisecondorX_reference.npz \
  > "$RES_INTERMEDIATE/${FILEID}.wisecondorx-b37_gender.txt"

fi

## Generate PDF-report
conda activate r_env

Rscript -e "rmarkdown::render('reportgeneration.Rmd', output_format='pdf_document', output_file='${RESULTSFOLDER}/${FILEID}/intermediateOutput/${FILEID}_report.pdf', output_dir='${RESULTSFOLDER}/${FILEID}/intermediateOutput', params=list(fileid='${FILEID}'), intermediates_dir='${RESULTSFOLDER}/${FILEID}/rmd_tmp')" "$FILEID"

mv *.log logs/.
