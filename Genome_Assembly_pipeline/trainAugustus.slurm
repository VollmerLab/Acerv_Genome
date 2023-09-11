#!/bin/bash

#SBATCH --job-name=Augustus
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#module load bedtools
module load singularity
MAKER=/work/vollmer/software/maker.sif

#outdir=/scratch/j.selwyn/tmp/maker_r1
outdir=${1}

mkdir -p ${outdir}/augustus
cd ${outdir}/augustus

# Convert files to be usable by augustus
ln -s ${outdir}/snap/export.ann ${outdir}/augustus/
ln -s ${outdir}/snap/export.dna ${outdir}/augustus/

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  zff2augustus_gbk.pl > ${outdir}/augustus/augustus.gbk

#Make testing/training split - training set will be 10% of genes or 100 loci whichever is smaller
n_loci=$(grep -c LOCUS augustus.gbk)
pct_loci=$(( 10*n_loci / 500 ))
test_loci=$(( pct_loci < 500 ? pct_loci : 500 ))
echo "Using ${test_loci} test loci out of a total of ${n_loci} loci"

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  randomSplit.pl ${outdir}/augustus/augustus.gbk ${test_loci}

#grep -c LOCUS augustus.gbk.train
#grep -c LOCUS augustus.gbk.test

#Create New Species to train for A cervicornis
export AUGUSTUS_CONFIG_PATH=${outdir}/../augustus_config/config

if [ ! -d "${outdir}/../augustus_config/config" ]
then
  mkdir -p ${outdir}/../augustus_config
  singularity exec --bind /work,/scratch,/tmp ${MAKER} \
    cp -r /usr/share/augustus/config ${outdir}/../augustus_config

  #Make New Augustus Species
  singularity exec --bind /work,/scratch,/tmp ${MAKER} \
    new_species.pl \
      --species=acropora_cervicornis
fi

#Train Augustus Part 1
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  etraining --species=acropora_cervicornis ${outdir}/augustus/augustus.gbk.train
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  augustus --species=acropora_cervicornis ${outdir}/augustus/augustus.gbk.test | tee ${outdir}/augustus/first_training.out


#Optimize Prediction Params
#default is 5 rounds and 8 kfold
#takes too long so drop rounds to 1 and repeat whole augustus process 5 times
#leave pstep as default 5. bump kfold to 15 to use more threads and maybe make
#each individual thread run faster
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  optimize_augustus.pl \
    --cpus=${SLURM_CPUS_PER_TASK} \
    --species=acropora_cervicornis \
     --rounds=1 \
     --kfold=15 \
     --pstep=5 \
    ${outdir}/augustus/augustus.gbk.train

#Second round of training
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  etraining --species=acropora_cervicornis ${outdir}/augustus/augustus.gbk.train
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  augustus --species=acropora_cervicornis ${outdir}/augustus/augustus.gbk.test | tee ${outdir}/augustus/second_training.out
