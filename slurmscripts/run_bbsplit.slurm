#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks=24
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export STRAIN='strain id'; sbatch --mail-user=$EMAIL -J "${STRAIN}-filter" scripts/run_bbsplit.slurm

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "FILTERING READS BASED ON BLOBTOOLS OUTPUT FOR $STRAIN"

ml conda-env/genomics-nw-py3.12.8
ml bbmap

ASMFILE="../../../assembly_12B1/ncbi/GCA_041296205.1_12B1_v1.1_genomic.fna"
OUTDIR="../../../hybrid_diversity/1_sequencing/3_filtered_subset_reads"
READ1="../../../hybrid_diversity/1_sequencing/2_filtered_reads/${STRAIN}_1.fq.gz"
READ2="../../../hybrid_diversity/1_sequencing/2_filtered_reads/${STRAIN}_2.fq.gz"




cd $OUTDIR

bbsplit.sh -Xmx16g nodisk=t in1=$READ1 in2=$READ2 ref=${STRAIN}_contamtigs.fa basename=%.fq.gz outu1=${STRAIN}_12B1_1.fq.gz outu2=${STRAIN}_12B1_2.fq.gz
reformat.sh samplereadstarget=30000000 in=${STRAIN}_12B1_1.fq in2=${STRAIN}_12B1_2.fq out=${STRAIN}_12B1_30m_1.fq out2=${STRAIN}_12B1_30m_2.fq
