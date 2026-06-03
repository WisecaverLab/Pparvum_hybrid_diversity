#!/bin/bash
#SBATCH --job-name=emapper
#SBATCH --account=jwisecav
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --time=14-00:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

#USAGE: sbatch --mail-user=$EMAIL scripts/run_eggnog.slurm

ml biocontainers
ml eggnog-mapper/2.1.7

emapper.py --cpu 32 -i GCA_041296205.1_12B1_v1.1_protein.faa 
