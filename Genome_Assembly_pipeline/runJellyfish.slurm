#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=jellyfish
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outdir=${1};shift 1 #Either DNA or RNA (changes how mapping is done due to splice sites)
fastq_files="$@"

mkdir -p ${outdir}
cd ${outdir}

MEMORY_ALOTTED=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 1" | bc)
echo "Memory allowed = ${MEMORY_ALOTTED}"
# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
sample_names=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/[\._][rR][12]\.*\..*//' | sort -u) )
extensions=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/^.*\([\._][rR][12]\)/\1/g' | sort -u) )
number_of_individuals=${#sample_names[@]}
number_of_extensions=${#extensions[@]}

sample_name=${sample_names[${SLURM_ARRAY_TASK_ID}]}

echo "Running Jellyfish on ${sample_name}"

module load anaconda3; source activate jellyfish
if (file ${indir}/${sample_name}${extensions[0]} | grep -q compressed ) ; then
  jellyfish \
    count \
      -C \
      -m 21 \
      -s ${MEMORY_ALOTTED}00000 \
      -t ${SLURM_CPUS_PER_TASK} \
      <(zcat ${indir}/${sample_name}${extensions[0]}) \
      <(zcat ${indir}/${sample_name}${extensions[1]}) \
      -o ${outdir}/${sample_name}.jf
  else
    jellyfish \
      count \
        -C \
        -m 21 \
        -s ${MEMORY_ALOTTED}00000 \
        -t ${SLURM_CPUS_PER_TASK} \
        ${indir}/${sample_name}${extensions[0]} \
        ${indir}/${sample_name}${extensions[1]} \
        -o ${outdir}/${sample_name}.jf
fi


jellyfish \
  histo \
    -t ${SLURM_CPUS_PER_TASK} \
    ${outdir}/${sample_name}.jf > ${outdir}/${sample_name}.histo

rm ${outdir}/${sample_name}.jf
