#!/bin/bash

#SBATCH --job-name=rds2h5ad
#SBATCH --partition=day
#SBATCH --time=24:00:00
#SBATCH --array=0
#SBATCH --cpus-per-task=20
#SBATCH --mem=384G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yoshiaki.yasumizu@yale.edu
#SBATCH --output=/home/yy693/pi_hafler/ASAP/script/slog/rds2h5ad-%A_%a.txt

# Cell Ranger output directories
cellranger_dirs=($(ls -d /home/yy693/pi_hafler/ASAP/output/cellranger/*))

# Get the directory for this array task
dir=${cellranger_dirs[$SLURM_ARRAY_TASK_ID]}
dirname=$(basename "$dir")

# Change to the working directory
cd /home/yy693/pi_hafler/ASAP

# make tmp dir in the scratch
mkdir -p /home/yy693/palmer_scratch/tmp/$SLURM_JOB_ID

R_LIBS_USER=/home/yy693/pi_hafler/R_singularity/R/bioconductor_docker_RELEASE_3_19-R-4.4.1
export R_LIBS_USER=${R_LIBS_USER}

# Run the script for this directory
apptainer exec  --contain \
  --bind /gpfs/gibbs/pi/hafler/yy693:/home/yy693/pi_hafler \
  --bind /home/yy693/palmer_scratch/tmp/$SLURM_JOB_ID:/tmp \
/vast/palmer/pi/hafler/yy693/R_singularity/bioconductor_docker_RELEASE_3_19-R-4.4.1.sif \
Rscript /home/yy693/pi_hafler/ASAP/script/241115_screfmapping_rds2h5ad.R "$dirname"