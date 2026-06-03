#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks=128
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export STRAIN='strain id'; export WIN='window size'; sbatch --mail-user=$EMAIL -J "${STRAIN}-${WIN}kb" scripts/run_bedtools.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "CALCULATING COVERAGE OF $WIN KB WINDOWS FOR $STRAIN"

CWD=$PWD
ml biocontainers
ml bedtools/2.31.0 

DIR="4_read_mapping/2_filtered_primary_reads/2_nuc_and_chloroplast"

cd $DIR

bedtools coverage -a ../../0_genome_reference/12B1_combined.windows_${WIN}kb.bed -b ${STRAIN}_to_12B1_sorted_filtered_mapq30.bam -mean > ${STRAIN}_to_12B1_sorted_filtered_mapq30.coverage_${WIN}kb.txt

cd $CWD


DIR="4_read_mapping/3_deduplicated/2_nuc_and_chloroplast"

cd $DIR

bedtools coverage -a ../../0_genome_reference/12B1_combined.windows_${WIN}kb.bed -b ${STRAIN}_markduplicates.bam -mean > ${STRAIN}_markduplicates.coverage_${WIN}kb.txt
