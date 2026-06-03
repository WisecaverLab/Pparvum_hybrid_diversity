#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export STRAIN='strain id'; sbatch --mail-user=$EMAIL -J "${STRAIN}-fastp" scripts/run_fastp_paired.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "RUNNING FASTP FOR $STRAIN"

#User defined parameters
INDIR='1_sequencing/2_filtered_reads'
OUTDIR='1_sequencing/2_filtered_reads'
THREADS=16
OPTIONS="--thread $THREADS --json ${OUTDIR}/${STRAIN}_fastp.json --html ${OUTDIR}/${STRAIN}_fastp.html fastp --length_required 150 --length_limit 151"

#module use /depot/jwisecav/apps/modulefiles       
#module load conda-env/read_preprocessing-py3.8.3 
#module load fastqc 

module load biocontainers
module load fastp/0.23.2

echo fastp -i ${INDIR}/${STRAIN}_clean_1.fq.gz -I ${INDIR}/${STRAIN}_clean_2.fq.gz -o ${OUTDIR}/${STRAIN}_cleanlenfilt_1.fq.gz -O ${OUTDIR}/${STRAIN}_cleanlenfilt_2.fq.gz $OPTIONS
time fastp -i ${INDIR}/${STRAIN}_clean_1.fq.gz -I ${INDIR}/${STRAIN}_clean_2.fq.gz -o ${OUTDIR}/${STRAIN}_cleanlenfilt_1.fq.gz -O ${OUTDIR}/${STRAIN}_cleanlenfilt_2.fq.gz $OPTIONS
