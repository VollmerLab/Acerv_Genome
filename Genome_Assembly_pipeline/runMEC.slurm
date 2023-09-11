#!/bin/bash

#SBATCH --job-name=MEC
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
module load samtools
module load bowtie/2.3.5.1

MEC=/work/vollmer/software/mec.sif
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
draft_genome=${2}; shift 2
polishing_reads=( "$@" )

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

mkdir -p ${outdir}
cd ${outdir}

# Map Reads to Genome
MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}


bowtie2-build ${draft_genome} contigs --threads ${SLURM_CPUS_PER_TASK}
bowtie2 -p ${SLURM_CPUS_PER_TASK} -x contigs -1 ${polishing_reads[0]} -2 ${polishing_reads[1]} -S contigs_short.sam

samtools view -Sb -@ ${SLURM_CPUS_PER_TASK} contigs_short.sam > contigs_short.bam
echo "Finished SAM -> BAM"
samtools sort -@ ${SLURM_CPUS_PER_TASK} -m ${MEMORY_PER_TASK}M -O BAM contigs_short.bam > contigs_short.sort.bam
echo "Finished BAM sorting"
samtools index -b -@ ${SLURM_CPUS_PER_TASK} contigs_short.sort.bam
echo "Finished BAM indexing"

# Break Poorly Covered Misassemblies
singularity exec --bind /work,/scratch,/tmp ${MEC} \
  python /opt/MEC/src/mec.py \
    -bam contigs_short.sort.bam \
    -i ${draft_genome} \
    -o ${outdir}/${genome_name}_mec.${draft_genome##*.} \
    -q 40 \
    -m 0 \
    -s 0 \
    -a 0.4 \
    -b 0.5 \
    -g 0.2

rm contigs_short.sam
rm contigs_short.bam
rm contigs_short.sort.bam*

# Spawn assembly check job
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_mec.${draft_genome##*.}
