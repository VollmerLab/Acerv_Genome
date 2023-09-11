#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=mapShort
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

draft_genome=${1}
forward_illumina=${2}
outDir=${3}

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

indir=${forward_illumina%/*}
final_portion=${forward_illumina##*/}
illuminaName=${final_portion%%-*}
tmp=${final_portion#*-}
preprocess_steps=${tmp%%.*}
extension=${tmp#*.}
new_extension=$(echo ${extension} | sed 's/R1/R2/g' | sed 's/r1/r2/g')

reverse_illumina=${indir}/${illuminaName}-${preprocess_steps}.${new_extension}

mkdir -p ${outDir}
cd ${outDir}

#ln -s ${draft_genome} ./
#draft_genome=${outDir}/${genome_name}.${draft_genome##*.}

#Rename Contigs to not have any whacky names
conda deactivate; module load anaconda3; source activate bbmap
rename.sh in=${draft_genome} out=${outDir}/${genome_name}.${draft_genome##*.} prefix=Acerv_scaffold
draft_genome=${outDir}/${genome_name}.${draft_genome##*.}
conda deactivate

TOTAL_MEMORY=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 100" | bc) #Amount of memory in GB
MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc) #Amount of memory in MB / N nodes
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

echo ${SLURM_MEM_PER_NODE}
echo ${TOTAL_MEMORY}
echo ${MEMORY_PER_TASK}


module load bwa
module load samtools

#Map reads to draft genome
bwa index ${draft_genome}
bwa mem \
  -t ${SLURM_CPUS_PER_TASK} \
  ${draft_genome} \
  ${forward_illumina} \
  ${reverse_illumina} > ${outDir}/${illuminaName}_aligned_${genome_name}.sam

samtools view -Sb -@ ${SLURM_CPUS_PER_TASK} ${outDir}/${illuminaName}_aligned_${genome_name}.sam > ${outDir}/${illuminaName}_aligned_${genome_name}.bam
echo "Finished SAM -> BAM"
# -m ${MEMORY_PER_TASK}M
samtools sort -@ ${SLURM_CPUS_PER_TASK} -O BAM ${outDir}/${illuminaName}_aligned_${genome_name}.bam > ${outDir}/${illuminaName}_aligned_${genome_name}_sorted.bam
echo "Finished BAM sorting"

samtools index -b -@ ${SLURM_CPUS_PER_TASK} ${outDir}/${illuminaName}_aligned_${genome_name}_sorted.bam
echo "Finished BAM indexing"

#Get Mapping Stats #not recent enough version to output as tsv
samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${outDir}/${illuminaName}_aligned_${genome_name}_sorted.bam > ${outDir}/mapping_stats

rm ${outDir}/${illuminaName}_aligned_${genome_name}.sam
rm ${outDir}/${illuminaName}_aligned_${genome_name}.bam
