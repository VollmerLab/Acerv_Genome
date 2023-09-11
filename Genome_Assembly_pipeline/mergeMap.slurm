#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=mergeMap
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load samtools/1.10

outdir=${1} #outdir must have "trinity" in it
ref_genome=${2} #assembled reference genome

mkdir -p ${outdir}
cd ${outdir}

genome_name=$(echo ${ref_genome##*/} | sed -e "s/.fasta$//" | sed -e "s/.fa$//")

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

samtools merge -@ ${SLURM_CPUS_PER_TASK} -O BAM ${outdir}/${genome_name}.bam ${outdir}/*.bam
echo "Finished BAM Merging!"

rm ${outdir}/*_sorted.bam*

samtools sort -@ ${SLURM_CPUS_PER_TASK} -m ${MEMORY_PER_TASK}M -O BAM ${outdir}/${genome_name}.bam > ${outdir}/${genome_name}_sorted.bam
echo "Finished BAM sorting!"

samtools index -b -@ ${SLURM_CPUS_PER_TASK} ${outdir}/${genome_name}_sorted.bam
echo "Finished BAM indexing!"

samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${outdir}/${genome_name}_sorted.bam > ${outdir}/rnaSeq_for_Trinity_mapping_stats
echo "Finished Everything!"

rm ${outdir}/${genome_name}.bam
