#!/bin/bash
#SBATCH --account kga
#SBATCH -c 1
#SBATCH --mem 1g
#SBATCH --time 1:00:00


########################
# Input argument
RUN_ID=$1
########################


source /home/$USER/miniforge3/etc/profile.d/conda.sh
conda activate nipt-snakemake-env

## Downsampling to 50M reads 
snakemake -s analysisScripts/scripts/snakemake-down-sample-fastq-files --config runID="${RUN_ID}" --unlock


## Alignment - seqFF - chrY 
snakemake -s analysisScripts/scripts/snakemake-nipt-v2.0 -p --config runID="${RUN_ID}" --unlock


## Wisecondor 
snakemake -s analysisScripts/scripts/snakemake-nipt-v2.0-wisecondor -p --config runID="${RUN_ID}" --unlock


## WisecondorX 
snakemake -s analysisScripts/scripts/snakemake-nipt-v2.0-wisecondorx -p --config runID="${RUN_ID}" --unlock


## Generate PDF-report 
snakemake -s analysisScripts/scripts/snakemake-create-pdf-report-v3.smk -p --config runID="${RUN_ID}" --unlock



conda deactivate 