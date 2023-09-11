#!/bin/bash

#SBATCH --job-name=gff2protein
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
AGAT=/work/vollmer/software/agat.sif

module load anaconda3; source activate pigz

outdir=${1}
indir=${2}

#outdir=/scratch/j.selwyn/phylogenetics/proteins
#indir=/scratch/j.selwyn/phylogenetics/raw_data

mkdir -p ${outdir}
cd ${outdir}

all_gff=$(ls ${indir}/*gff*)
all_fasta=( $(echo ${all_gff} | sed 's/gff/fasta/g') )
all_gff=( ${all_gff} )

for (( i=0; i<=${#all_gff[@]}; i++ ))
do
  gff_file=${all_gff[i]}
  fasta_file=${all_fasta[i]}

  pigz -p ${SLURM_CPUS_PER_TASK} -d ${gff_file}
  pigz -p ${SLURM_CPUS_PER_TASK} -d ${fasta_file}

  gff_file=${gff_file%*.gz}
  fasta_file=${fasta_file%*.gz}
  outname=${fasta_file##*/}

  echo "###########"
  echo Extracting ${outname}
  echo "###########"

  singularity exec --bind /work,/scratch,/tmp ${AGAT} \
    agat_sp_extract_sequences.pl \
      --gff ${gff_file} \
      --fasta ${fasta_file} \
      --protein \
      -o ${outdir}/${outname} \
      --clean_final_stop \
      --clean_internal_stop

    pigz -p ${SLURM_CPUS_PER_TASK} ${gff_file}
    pigz -p ${SLURM_CPUS_PER_TASK} ${fasta_file}
done
