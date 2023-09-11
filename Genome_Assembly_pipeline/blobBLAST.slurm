#!/bin/bash

#SBATCH --job-name=blobBlast
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

#Blast Genome against database
module load anaconda3; source activate diamond
diamond \
  blastx \
    --threads ${SLURM_CPUS_PER_TASK} \
    --verbose \
    --query ${draft_genome} \
    --db ${protein_database} \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive \
    --max-target-seqs 5 \
    --evalue 1e-25 > ${outdir}/diamondBlastx.out

#Add blast hits to blob
singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
  blobtools \
    add \
      --hits ${outdir}/diamondBlastx.out \
      --taxrule bestsumorder \
      --taxdump ${taxdump_directory} \
      --threads ${SLURM_CPUS_PER_TASK} \
      ${outdir}/${genome_name}_blob
