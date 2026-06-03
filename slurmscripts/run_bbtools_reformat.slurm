#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export READS=8000000; export STRAIN='strain id'; sbatch --mail-user=$EMAIL -J "${STRAIN}-bbtools" scripts/run_bbtools_reformat.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "RUNNING BBTOOLS REFORMAT FOR $STRAIN"

#User defined parameters
INDIR='1_sequencing/2_filtered_reads'
OUTDIR='1_sequencing/3_filtered_subset_reads'

module load biocontainers
module load bbmap

reformat.sh reads=$READS sampleseed=12345 \
in1=${INDIR}/${STRAIN}_cleanlenfilt_1.fq.gz \
in2=${INDIR}/${STRAIN}_cleanlenfilt_2.fq.gz \
out1=${OUTDIR}/${STRAIN}_subsample${READS}_1.fq.gz \
out2=${OUTDIR}/${STRAIN}_subsample${READS}_2.fq.gz   

