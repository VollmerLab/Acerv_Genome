#!/bin/bash

#SBATCH --job-name=GeneMark
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
MAKER=/work/vollmer/software/maker.sif

#outdir=/scratch/j.selwyn/tmp/maker_r1
#genome=/scratch/j.selwyn/tmp/k2_flyeHQ_polish_1000000.fa

outdir=${1}
genome=${2}

mkdir -p ${outdir}/genemark
cd ${outdir}/genemark

singularity exec --bind /work,/scratch,/tmp ${MAKER} cp /usr/local/Genemark/gm_key ~/.gm_key

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  gmes_petap.pl \
    --ES \
    --v \
    --cores ${SLURM_CPUS_PER_TASK} \
    --sequence ${genome}


#output need to go
cp ${outdir}/genemark/gmhmm.mod ${outdir}/gmhmm.mod
