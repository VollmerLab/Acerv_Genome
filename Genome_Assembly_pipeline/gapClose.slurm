#!/bin/bash

#SBATCH --job-name=gapClose
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G

outdir=${1}
draft_genome=${2}
long_polishing_reads=${3}
#short_read_only_assembly=${4} - use short read assembly for gap filling too??

module load minimap2

scriptDir=/work/vollmer/software/jds_scripts

USABLE_MEMORY=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 1000" | bc)
echo "Memory call = ${USABLE_MEMORY}"

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

mkdir -p ${outdir}
cd ${outdir}

#Reformat input reads from fastq to fasta
module load anaconda3; source activate seqtk
seqtk seq -a ${long_polishing_reads} > ${outdir}/nanopore_reads.fasta

# Add in error correction with short reads?? Probably unnecessary due to already using racon & low error rate with new chemistry
conda deactivate; module load anaconda3; source activate tgsgapcloser
tgsgapcloser \
  --scaff ${draft_genome} \
  --reads ${outdir}/nanopore_reads.fasta \
  --output ${outdir}/gapClosed_assembly \
  --thread ${SLURM_CPUS_PER_TASK} \
  --tgstype ONT \
  --ne
mv ${outdir}/gapClosed_assembly.scaff_seq* ${outdir}/${genome_name}_gapClose.${draft_genome##*.}

# Spawn assembly check job
conda deactivate
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_gapClose.${draft_genome##*.}
