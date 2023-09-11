#!/bin/bash

#SBATCH --job-name=blobStart
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
BLOBTOOLS=/work/vollmer/software/blobtools2.sif

outdir=${1}
draft_genome=${2}

taxid=6130 #https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6130
protein_database=/work/vollmer/Databases/blobtools/uniprot/reference_proteomes.dmnd
taxdump_directory=/work/vollmer/Databases/blobtools/taxdump

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

mkdir -p ${outdir}
cd ${outdir}

#Initial BlobTools setup
singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
  blobtools \
    create \
      --fasta ${draft_genome} \
      --taxid ${taxid} \
      --taxdump ${taxdump_directory} \
      --threads ${SLURM_CPUS_PER_TASK} \
      ${outdir}/${genome_name}_blob
