#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=gatherPilon
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G

#Script to assemble polish a draft assembly using illumina sequences

#forward_illumina=/scratch/j.selwyn/genome_assemblies/k2_hybrid_28-July-2022/repair/k2_illumina-fp-fqscrn-repr.R1.fastq.gz

draft_genome=${1}
forward_illumina=${2}
outDir=${3}
scriptDir=/work/vollmer/software/jds_scripts

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

indir=${forward_illumina%/*}
final_portion=${forward_illumina##*/}
illuminaName=${final_portion%%-*}
tmp=${final_portion#*-}
preprocess_steps=${tmp%%.*}
extension=${tmp#*.}
new_extension=$(echo ${extension} | sed 's/R1/R2/g' | sed 's/r1/r2/g')

reverse_illumina=${indir}/${illuminaName}-${preprocess_steps}.${new_extension}

TOTAL_MEMORY=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 1000" | bc)
MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

echo ${TOTAL_MEMORY}
echo ${MEMORY_PER_TASK}

mkdir -p ${outDir}
cd ${outDir}

ln -s ${draft_genome} ./
draft_genome=${outDir}/${genome_name}.${draft_genome##*.}

#Find all pilon out and join together
cat ${outDir}/${genome_name}_*_pilon.fasta > ${outDir}/${genome_name}_pilon.fasta

sbatch \
  --output=${outDir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outDir}/${genome_name}_pilon.fasta
