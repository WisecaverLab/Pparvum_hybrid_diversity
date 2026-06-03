#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks=12
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export STRAIN='strain id'; sbatch --mail-user=$EMAIL -J "${STRAIN}-${WIN}kb" scripts/run_samtools.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "CALCULATING COVERAGE OF $WIN KB WINDOWS FOR $STRAIN"

CWD=$PWD
ml biocontainers
ml samtools/1.22.1 

DIR="4_read_mapping/3_deduplicated/2_nuc_and_chloroplast"

cd $DIR

samtools index ${STRAIN}_markduplicates.bam

samtools view -b ${STRAIN}_markduplicates.bam $(cat nuclear_chroms.txt | tr '\n' ' ') > ${STRAIN}_markduplicates_nuconly.bam

samtools index ${STRAIN}_markduplicates_nuconly.bam
