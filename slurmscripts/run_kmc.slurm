#!/bin/bash
#SBATCH --account=jwisecav
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00     
#SBATCH --output=slurm-out/%x.o%j
#SBATCH --error=slurm-out/%x.e%j

# usage: export READS=30; export STRAIN='strain id'; sbatch --mail-user=$EMAIL -J "${STRAIN}-kmc" scripts/run_kmc.slurm 

[  -z "$STRAIN" ] && echo "STRAIN NOT SPECIFIED!" && exit || echo "CREATING KMC HISTOGRAMS FOR $STRAIN"

INDIR="1_sequencing/3_filtered_subset_reads"
OUTDIR="3_kmer_frequencies/"
TIME=`date +'%s'`
TMPDIR="/tmp/hetplot_${STRAIN}_${TIME}"
KMER=21
MAXKMER=10000
THREADS=30

echo "Creating working directories"
if [ ! -d "$OUTDIR" ]; then mkdir -p $OUTDIR; fi

# Load modules
module load biocontainers       
module load kmc/3.2.1

# Create kmc histogram of read depth
mkdir $TMPDIR 
ls ${INDIR}/${STRAIN}_12B1_${READS}m_*.fq.gz > ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}.files

echo -e "kmc -k21 -t16 -m60 -ci1 -cs${MAXKMER} @${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}.files ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}_reads $TMPDIR\n"
kmc -k21 -t16 -m60 -ci1 -cs${MAXKMER} @${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}.files ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}_reads $TMPDIR

echo -e "kmc_tools transform ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}_reads histogram ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}.histo -cx${MAXKMER}\n"
kmc_tools transform ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}_reads histogram ${OUTDIR}/${STRAIN}_12B1_${READS}m_k${KMER}.histo -cx${MAXKMER}

