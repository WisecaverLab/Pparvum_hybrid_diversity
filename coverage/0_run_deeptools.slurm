#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks=2
#SBATCH --time=2:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export STRAIN='strain id'; export BIN=10000; sbatch --mail-user=$EMAIL -J "${STRAIN}-cov${BIN}" scripts/run_deeptools.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "CALCULATING COVERAGE OF $BIN WINDOWS FOR $STRAIN"

CWD=$PWD
ml biocontainers
ml deeptools/3.5.5

DIR="4_read_mapping/3_deduplicated/2_nuc_and_chloroplast"
KBP=$((BIN / 1000))
ASMLEN=93538723

cd $DIR

echo bamCoverage -b ${STRAIN}_markduplicates_nuconly.bam -o ${STRAIN}_markduplicates_nuconly.coverage_${KBP}kb.bw --binSize $BIN --normalizeUsing RPGC --effectiveGenomeSize $ASMLEN
bamCoverage -b ${STRAIN}_markduplicates_nuconly.bam -o ${STRAIN}_markduplicates_nuconly.coverage_${KBP}kb.bw --binSize $BIN --normalizeUsing RPGC --effectiveGenomeSize $ASMLEN

