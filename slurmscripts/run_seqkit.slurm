#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: sbatch --mail-user=$EMAIL -J "seqkit" scripts/run_seqkit.slurm

ml biocontainers seqkit

cd 1_sequencing/2_filtered_reads/
seqkit stats --threads 30 *.fq.gz

