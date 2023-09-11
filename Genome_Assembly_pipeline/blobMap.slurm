#!/bin/bash

#SBATCH --job-name=blobMap
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load minimap2/2.17
module load samtools
module load singularity
BLOBTOOLS=/work/vollmer/software/blobtools2.sif

outdir=${1}
draft_genome=${2}
long_polishing_reads=${3}; shift 3
short_polishing_reads=( "$@" )

taxid=6130 #https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6130
protein_database=/work/vollmer/Databases/blobtools/uniprot/reference_proteomes.dmnd
taxdump_directory=/work/vollmer/Databases/blobtools/taxdump

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

mkdir -p ${outdir}
cd ${outdir}

#Map Long Reads to Genome
minimap2 -ax sr \
         -t ${SLURM_CPUS_PER_TASK} \
         ${draft_genome} \
         ${long_polishing_reads} \
| samtools sort -@${SLURM_CPUS_PER_TASK} -O BAM -o ${outdir}/long_reads.bam -
echo "Finished sorting long reads"

#Map Short Reads to Genome
minimap2 -ax sr \
         -t ${SLURM_CPUS_PER_TASK} \
         ${draft_genome} \
         ${short_polishing_reads[@]} \
| samtools sort -@${SLURM_CPUS_PER_TASK} -O BAM -o ${outdir}/short_reads.bam -
echo "Finished sorting short reads"

#Merge mapped reads
samtools merge -@ ${SLURM_CPUS_PER_TASK} ${outdir}/merged_reads.bam ${outdir}/long_reads.bam ${outdir}/short_reads.bam

rm ${outdir}/short_reads.bam
rm ${outdir}/long_reads.bam

#Add mapped reads to blob
singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
  blobtools \
    add \
      --cov ${outdir}/merged_reads.bam \
      --threads ${SLURM_CPUS_PER_TASK} \
      ${outdir}/${genome_name}_blob
