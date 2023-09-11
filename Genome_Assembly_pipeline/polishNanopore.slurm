#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=polish_nanopore
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
module load minimap2/2.17

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

minimap2 \
  -a \
  -x map-ont \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam \
  ${draft_genome} \
  ${polishing_reads}

#Polish with Racon
module load anaconda3; source activate racon

racon \
  -u \
  -m 3 \
  -x -5 \
  -g -4 \
  -w 500 \
  -t ${SLURM_CPUS_PER_TASK} \
  ${polishing_reads} \
  ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam \
  ${draft_genome} > ${outDIR}/${genome_name}_racon.fasta

rm ${outDIR}/${genome_name}_aligned_${nanopore_name}.sam

# Spawn assembly check job
conda deactivate
scriptDir=/work/vollmer/software/jds_scripts
sbatch \
  --output=${outDIR}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outDIR}/${genome_name}_racon.${draft_genome##*.}
