#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:v100-sxm2:1
#SBATCH --time=08:00:00
#SBATCH --job-name=guppy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=180G

#Script used to take folder of fast5 files from ONT Mk1C sequencer and perform basecalling. Additionally produces several graphical and summary statistical outputs about the raw basecalled reads

outDir=${1}; shift 1
sequencingPaths=( "$@" ) #path containing all input files within - will recursively search to find all fast5 files

sequencingPath=${sequencingPaths[${SLURM_ARRAY_TASK_ID}]}

outName=${sequencingPath##*/}

mkdir -p ${outDir}/${SLURM_ARRAY_TASK_ID}

module load guppy/6.1.7

#Edit below with other settings desired for basecalling
guppy_basecaller \
  -i ${sequencingPath} \
  -r \
  -s ${outDir}/${SLURM_ARRAY_TASK_ID} \
  --flowcell FLO-MIN112 \
  --kit SQK-LSK112 \
  -x cuda:all \
  --verbose_logs

#Concatenate all FASTQ files into one
find ${outDir}/${SLURM_ARRAY_TASK_ID} -type f -name '*.fastq' -exec cat {} + > ${outDir}/${outName}.fastq

#Run some checks and basic stats
#module load anaconda3; source activate fastq_utils
#fastq_info ${outDir}/${outName}.fastq

#Currently not working (17-June-2022)
#conda deactivate; module load anaconda3; source activate nanostat
#NanoStat \
#  --fastq ${outDir}/${outName}.fastq \
#  -t ${SLURM_CPUS_PER_TASK} \
#  --name NanoStat_Summary_report \
#  --outdir ${outDir}/nanostat/ \
#  --tsv

#Currently only partly working (17-June-2022)
#conda deactivate; module load anaconda3; source activate nanoplot
#NanoPlot \
#  -t ${SLURM_CPUS_PER_TASK} \
#  --fastq ${outDir}/${outName}.fastq \
#  --plots hex dot \
#  --title ${outName} \
#  --outdir ${outDir}/nanoplot/

#Compress assembly
module load anaconda3; source activate pigz
#gzip ${outDir}/${outName}.fastq
pigz -p ${SLURM_CPUS_PER_TASK} "${outDir}/${outName}.fastq"

#To run:
#sbatch scripts/nanopore/basecallGuppy.surm \
#   /work/vollmer/Nanopore_Data/k2cerv_1_16-June-2022/k2cerv
#   /work/vollmer/Nanopore_Data/k2cerv_1_16-June-2022/baseCalls
#   k2cerv
