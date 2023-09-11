#!/bin/bash

#SBATCH --job-name=SNAP
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
MAKER=/work/vollmer/software/maker.sif

#outdir=/scratch/j.selwyn/tmp/maker_r3
#gff_in=/scratch/j.selwyn/tmp/maker_r3/r3.gff

outdir=${1}
gff_in=${2}

mkdir -p ${outdir}/snap/params
cd ${outdir}/snap

#### Get some stats about the previous round ####
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  maker2zff \
    -n \
    ${gff_in}

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  fathom \
    ${outdir}/snap/genome.ann \
    ${outdir}/snap/genome.dna \
    -gene-stats > ${outdir}/gene-stats.log 2>&1


#### Pull out high confidence genes, validate, and train model
# export 'confident' gene models from MAKER and rename to something meaningful
#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  maker2zff \
#    -x 0.25 \
#    -l 50 \
#    ${gff_in}

# have default filtering settings following https://reslp.github.io/blog/My-MAKER-Pipeline/
# this will overwrite the *ann and *dna files which contain all proteins with only those which meet QC controls
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  maker2zff \
    -n \
    ${gff_in}

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  fathom \
    ${outdir}/snap/genome.ann \
    ${outdir}/snap/genome.dna \
    -validate > ${outdir}/validate.log 2>&1

#Find & remove errors in gene models
cat ${outdir}/validate.log | grep "error" > ${outdir}/snap_model_errors.txt
n_errors=$(grep -c "MODEL" ${outdir}/snap_model_errors.txt)
error_models=$(grep -Po "MODEL[0-9]+" ${outdir}/snap_model_errors.txt)

echo "Removing ${n_errors} models with errors"
for model in ${error_models}
do
  echo "Removing model ${model}"
  grep -vwE ${model} ${outdir}/snap/genome.ann > ${outdir}/snap/tmp_file
  mv ${outdir}/snap/tmp_file ${outdir}/snap/genome.ann
done

# collect the training sequences and annotations, plus 1000 surrounding bp for training
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  fathom \
    ${outdir}/snap/genome.ann \
    ${outdir}/snap/genome.dna \
    -categorize 1000 > categorize.log 2>&1

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  fathom \
    ${outdir}/snap/uni.ann \
    ${outdir}/snap/uni.dna \
    -export 1000 \
    -plus > uni-plus.log 2>&1

# create the training parameters
cd ${outdir}/snap/params
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  forge \
    ${outdir}/snap/export.ann \
    ${outdir}/snap/export.dna > ${outdir}/snap/forge.log 2>&1
cd ${outdir}

# assembly the HMM
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  hmm-assembler.pl \
    genome \
    snap/params > ${outdir}/snap_modelTrain.hmm
