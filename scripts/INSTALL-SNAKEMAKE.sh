#!/bin/bash


source "/home/$USER/miniforge3/etc/profile.d/conda.sh"

conda create --name snakemake_env
conda activate snakemake_env
conda install bioconda::snakemake -y
conda deactivate



