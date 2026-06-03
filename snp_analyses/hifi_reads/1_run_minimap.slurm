#!/bin/bash
#SBATCH --job-name minimap
#SBATCH --partition=wisecaver
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time 12:00:00
#SBATCH --output=slurm-out/%x.o%j             
#SBATCH --error=slurm-out/%x.e%j

# usage: sbatch scripts/snp_analyses/hifi_reads/run_minimap.slurm

module load minimap2/2.26
module load samtools/1.21

THREADS=16
GENOME="4_read_mapping/0_genome_reference/GCA_041296205.1_12B1_v1.1_genomic.fna"
READS="../../data/cngb/UTEX2797_hifi.fq"
OPTIONS="-ax map-hifi --secondary=no"
KMER=13

OUTFILE="4_read_mapping/4_hifi_reads/map-hifi/UTEX2797_hifi_v12B1_K${KMER}"

minimap2 -t $THREADS $OPTIONS -k$KMER $GENOME $READS > ${OUTFILE}.sam
samtools view -@ $THREADS -bS -o ${OUTFILE}.bam ${OUTFILE}.sam
samtools sort -@ $THREADS -o ${OUTFILE}.sorted.bam ${OUTFILE}.bam
samtools index -@ $THREADS ${OUTFILE}.sorted.bam

samtools view -b -F 4 -F 2048 -q 20 ${OUTFILE}.sorted.bam > ${OUTFILE}_filtered_q20.bam
samtools index ${OUTFILE}_filtered_q20.bam
