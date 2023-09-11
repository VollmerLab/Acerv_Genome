#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=flye
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=180G

#Script to assemble nanopore fastq file into a draft genome using the Flye assembly pipeline https://github.com/fenderglass/Flye

nanopore_file=${1}
outDir=${2}
genomeSize=${3} #can use k/m/g suffix e.g. 500m
genome_name=${4}

mkdir -p ${outDir}
cd ${outDir}

module load anaconda3
source activate flye

#Modify code after this for all the Flye parameters
flye \
  --nano-hq ${nanopore_file} \
  --read-error 0.03 \
  -g ${genomeSize} \
  -o ${outDir} \
  -t ${SLURM_CPUS_PER_TASK}

mv ${outDir}/assembly.fasta ${outDir}/${genome_name}_flyeHQ.fasta

# Spawn assembly check job
scriptDir=/work/vollmer/software/jds_scripts
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outDir}/${genome_name}_flyeHQ.fasta

#To run:

#sbatch scripts/nanopore/assembleShasta.slurm \
#    ~/test_assembly/drosphilla/SRR6702603_1.fastq \
#    ~/test_assembly/flye/drosphilla
#    180m
