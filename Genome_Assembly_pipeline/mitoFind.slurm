#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=mitoFind
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
CONTAINER=/work/vollmer/software/mitoz.sif

array_id=${SLURM_ARRAY_TASK_ID}
outdir=${1}; shift 1
fastq_files="$@"

USABLE_MEMORY=$(echo "scale=0 ; 0.8 * ${SLURM_MEM_PER_NODE} / 1000" | bc)

# Sort out post-array & choose sample to run in array
tmp=( ${fastq_files} )
indir=${tmp%/*}
sample_names=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/[\._][rR][12]\.*\..*//' | sort -u) )
extensions=( $(ls ${fastq_files} | sed -e 's!.*/!!' | sed -e 's/^.*\([\._][rR][12]\)/\1/g' | sort -u) )
number_of_individuals=${#sample_names[@]}
number_of_extensions=${#extensions[@]}

sample_name=${sample_names[${array_id}]}

echo "Assembling ${sample_name} Mitochondrial Genome"

mkdir -p ${outdir}/${sample_name}
cd ${outdir}/${sample_name}

#Trick it into using Cnidarian proteins
#https://github.com/linzhi2013/MitoZ/issues/181
mkdir -p ${outdir}/${sample_name}/custom_profile
singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
  cp -a \
    /app/anaconda/lib/python3.9/site-packages/mitoz-3.5-py3.9.egg/mitoz/profiles \
    ${outdir}/${sample_name}/custom_profile
 mv ${outdir}/${sample_name}/custom_profile/profiles/MT_database/Arthropoda_CDS_protein.fa \
  ${outdir}/${sample_name}/custom_profile/profiles/MT_database/bak.Arthropoda_CDS_protein.fa
 ln -s ${outdir}/${sample_name}/custom_profile/profiles/MT_database/Cnidaria_CDS_protein.fa \
  ${outdir}/${sample_name}/custom_profile/profiles/MT_database/Arthropoda_CDS_protein.fa


#mitoZ
singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
  mitoz \
    all \
      --profiles_dir ${outdir}/${sample_name}/custom_profile/profiles \
      --clade Arthropoda \
      --outprefix ${sample_name} \
      --workdir ${outdir}/${sample_name} \
      --thread_number ${SLURM_CPUS_PER_TASK} \
      --species_name acropora_${sample_name} \
      --skip_filter \
      --fq1 ${indir}/${sample_name}${extensions[0]} \
      --fq2 ${indir}/${sample_name}${extensions[1]} \
      --assembler megahit \
      --memory ${USABLE_MEMORY} \
      --requiring_taxa Cnidaria

#Clean-up
