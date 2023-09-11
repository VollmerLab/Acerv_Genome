#!/bin/bash

#SBATCH --job-name=RepeatModeler
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G

#specifically warned against breaking up genome and running on multiple nodes

module load singularity
MAKER=/work/vollmer/software/maker.sif

#genome=/work/vollmer/k2_nanopore_genome/k2_flyeHQ_polish.fa
#outdir=/scratch/j.selwyn/tmp/repeatModelerFull

outdir=${1}
genome=${2}


mkdir -p ${outdir}
cd ${outdir}

genome_name=$(echo ${genome##*/} | sed -e "s/\.f.*a$//")

#module load anaconda3; source activate bioawk
#bioawk -v LENGTH="${minLength_genome}" -c fastx '{ if(length($seq) > LENGTH) { print ">"$name; print $seq }}' \
#  ${genome} > ${outdir}/${genome_name}_${minLength_genome}.fa
#conda deactivate

#Repeat Modeller
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  BuildDatabase \
    -name ${genome_name} \
    -engine ncbi \
    ${genome}

#https://github.com/Dfam-consortium/RepeatModeler/issues/158
export BLAST_USAGE_REPORT=FALSE
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  RepeatModeler \
    -pa $(( SLURM_CPUS_PER_TASK / 4 )) \
    -engine ncbi \
    -database ${genome_name} 2>&1 | tee ${outdir}/repeatmodeler.log
