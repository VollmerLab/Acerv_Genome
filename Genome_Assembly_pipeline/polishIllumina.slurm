#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=polishIllumina
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#Script to assemble polish a draft assembly using illumina sequences

#forward_illumina=/scratch/j.selwyn/genome_assemblies/k2_hybrid_28-July-2022/repair/k2_illumina-fp-fqscrn-repr.R1.fastq.gz

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

TOTAL_MEMORY=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 1000" | bc)
MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

echo ${TOTAL_MEMORY}
echo ${MEMORY_PER_TASK}

mkdir -p ${outDir}
cd ${outDir}

ln -s ${draft_genome} ./
draft_genome=${outDir}/${genome_name}.${draft_genome##*.}

#Split into which contig(s) to polish given array ID
split_jobs () {
  #https://www.geeksforgeeks.org/split-the-number-into-n-parts-such-that-difference-between-the-smallest-and-the-largest-part-is-minimum/
  X=$1
  N=$2

  out=""

  if [ $((X % N)) -eq 0 ]
    then
      for ((i=1; i<=N; i++))
      do
        out="${out} $((X / N))"
      done

    else
      zp=$((N - (X % N)))
      pp=$((X / N))

      for ((i=1; i<=N; i++))
      do
        if [ ${i} -gt ${zp} ]
        then
          out="${out} $((pp + 1))"
        else
          out="${out} ${pp}"
        fi
      done
    fi
  echo "${out}"
}

#ARRAY_ID=0
array_size=100
ARRAY_ID=${SLURM_ARRAY_TASK_ID}

all_contigs=( $(grep -o '^>.*' ${draft_genome} | sed -e 's/>//g') )
contig_bins=( $(split_jobs ${#all_contigs[@]} ${array_size}) )
starting_index=$(echo "0 ${contig_bins[@]:0:$((ARRAY_ID))} 0" | sed -e 's/  */+/g' | bc)
contig_names=( ${all_contigs[@]:${starting_index}:${contig_bins[${ARRAY_ID}]}} )
printf -v contig_names '%s,' "${contig_names[@]}"
contig_names="${contig_names%,}"

#Polish with Pilon
module load anaconda3; source activate pilon
pilon -Xmx${TOTAL_MEMORY}G \
  --genome ${draft_genome} \
  --targets ${contig_names} \
  --frags ${outDir}/${illuminaName}_aligned_${genome_name}_sorted.bam \
  --output ${genome_name}_${ARRAY_ID}_pilon \
  --outdir ${outDir} \
  --changes \
  --threads ${SLURM_CPUS_PER_TASK}

#rm ${outDir}/${illuminaName}_aligned_${genome_name}_sorted.bam*

# Spawn assembly check job
#scriptDir=/work/vollmer/software/jds_scripts
#sbatch \
#  --output=${outdir}/assembly_assessment_%j.output \
#  ${scriptDir}/assessAssembly.slurm \
#    ${outdir}/${genome_name}_pilon.${draft_genome##*.}
