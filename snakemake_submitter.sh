#!/bin/bash
#
#SBATCH --job-name=<job_name>
#SBATCH --ntasks=1   
#SBATCH --partition=compute
#SBATCH --time=48:00:00
#SBATCH --mem=2gb
#SBATCH --output=<path/to/work/directory>/log/<job_name>_snakefile_.%j.txt
#SBATCH --error=<path/to/work/directory>/log/<job_name>_snakefile_.%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=<user@clemson.edu>

cd <path/to/working/directory>
#mkdir -p ./{log,logs_slurm}

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
conda activate snakemake

#--dag | display | dot
#-p -n \

snakemake \
-s Snakefile \
--profile slurm \
--latency-wait 150 \
-p \
