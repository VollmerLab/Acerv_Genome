#!/bin/bash

#SBATCH --job-name=faBLAST
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#proteins=/scratch/j.selwyn/tmp/r3.all.maker.proteins.fasta
#array_size=100
#ARRAY_ID=0

module load anaconda3

outdir=${1}
proteins=${2}
array_size=${3}
ARRAY_ID=${SLURM_ARRAY_TASK_ID}

mkdir -p ${outdir}
cd ${outdir}

#### figure out which subset of proteins to use here
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

protein_names=$(grep -o '^>.*' ${proteins} | sed -e 's/>//g' | sed -e 's/ /@/g')
protein_names=( ${protein_names} )
n_prot=${#protein_names[@]}
array_size=$(( n_prot < array_size ? n_prot : array_size ))


protein_bins=( $(split_jobs ${#protein_names[@]} ${array_size}) )
starting_index=$(echo "0 ${protein_bins[@]:0:$((ARRAY_ID))} 0" | sed -e 's/  */+/g' | bc)
protein_names_to_use=( ${protein_names[@]:${starting_index}:${protein_bins[${ARRAY_ID}]}} )
echo "KOFAMscan Proteins: ${protein_names_to_use[@]}"

tmpDir=${outdir}/kofamscan_array/${ARRAY_ID}
mkdir -p ${tmpDir}
cd ${tmpDir}

#Subset the proteome to the relevant portions
module load anaconda3; source activate seqtk
#printf '%s\n' "${genome_contig_names[@]}" | seqtk subseq ${genome} - > ${tmpDir}/genome_contigs.fasta
touch ${tmpDir}/protein_chunks.fasta
for contig in ${protein_names_to_use[@]}; do

  contig_tweak=$(echo ${contig} | sed -e 's/@/ /g')
  echo ${contig_tweak}
  printf "${contig_tweak}\n" | seqtk subseq ${proteins} - >> ${tmpDir}/protein_chunks.fasta
done
conda deactivate

#Run annotations with subset of proteins
module load anaconda3; source activate kofamscan
exec_annotation \
  -o ${tmpDir}/output.kofamscan \
  --profile /work/vollmer/Databases/kegg/profiles \
  --ko-list /work/vollmer/Databases/kegg/ko_list \
  --cpu ${SLURM_CPUS_PER_TASK} \
  -f mapper \
  ${tmpDir}/protein_chunks.fasta
