#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=polish_medaka
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

draft_genome=${1}
polishing_reads=${2}
outDIR=${3}

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\\.f.*a$//")
nanopore_name=$(echo ${polishing_reads##*/} | sed -e "s/.fastq.gz//")

mkdir -p ${outDIR}
cd ${outDIR}

ln -s ${draft_genome} ./
draft_genome=${outDIR}/${genome_name}.${draft_genome##*.}

#map polishing reads to assembly
module load minimap2/2.17; module load samtools

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
#MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * 100000 / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

minimap2 \
  -a \
  -x map-ont \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam \
  ${draft_genome} \
  ${polishing_reads}

samtools view -Sb -@ ${SLURM_CPUS_PER_TASK} ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam > ${outDIR}/${genome_name}_aligned_${nanopore_name}.bam
echo "Finished SAM -> BAM"

samtools sort -@ ${SLURM_CPUS_PER_TASK} -m ${MEMORY_PER_TASK}M -O BAM ${outDIR}/${genome_name}_aligned_${nanopore_name}.bam > ${outDIR}/${genome_name}_aligned_${nanopore_name}_sorted.bam
echo "Finished BAM sorting"

samtools index -b -@ ${SLURM_CPUS_PER_TASK} ${outDIR}/${genome_name}_aligned_${nanopore_name}_sorted.bam
echo "Finished BAM indexing"

#Polish with MEDAKA
module load singularity
CONTAINER=/work/vollmer/software/medaka.sif
#singularity exec --bind /work,/scratch,/tmp ${CONTAINER} medaka --help

singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
  medaka \
    consensus \
      ${outDIR}/${genome_name}_aligned_${nanopore_name}_sorted.bam \
      ${outDIR}/${genome_name}.hdf \
      --model r941_min_hac_g507 \
      --threads ${SLURM_CPUS_PER_TASK}

#
singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
  medaka \
    stitch \
      --threads ${SLURM_CPUS_PER_TASK} \
      ${outDIR}/${genome_name}.hdf \
      ${draft_genome} \
      ${outDIR}/${genome_name}_medaka.fasta

#Cleanup
rm ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam
rm ${outDIR}/${genome_name}_aligned_${nanopore_name}.bam
rm ${outDIR}/${genome_name}_aligned_${nanopore_name}_sorted.bam*
rm ${outDIR}/${genome_name}.hdf


# Spawn assembly check job
scriptDir=/work/vollmer/software/jds_scripts
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_medake.${draft_genome##*.}
