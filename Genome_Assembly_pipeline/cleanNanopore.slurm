#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=clean_nanopore
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G


#Script to trim and filter nanopore reads using porechop & nanofilt

in_fastq=${1} #raw nanopore reads coming in
outDIR=${2} #directory to write output to - outputs will be automatically labelled

#Set up variable names
zipped_name=${in_fastq##*/}
extension=${zipped_name#*.}
base_name=${zipped_name%%.*}

porechop_out=${outDIR}/${base_name}_pc.${extension}
nanofilt_assemble=${outDIR}/${base_name}_pc_nA.${extension} #longer reads of lower quality for assembly
nanofilt_polish=${outDIR}/${base_name}_pc_nP.${extension} #shorter reads of higher quality for polishing

echo Outputs will be:
echo ${nanofilt_assemble}
echo ${nanofilt_polish}

mkdir -p ${outDIR}
cd ${outDIR}
#Run porechop first
module load anaconda3; source activate porechop

porechop \
  -i ${in_fastq} \
  -o ${porechop_out} \
  --threads ${SLURM_CPUS_PER_TASK}

#Run Trimfilt
conda deactivate; module load anaconda3; source activate nanofilt

#Taken from R. Gatins https://github.com/remygatins/Holacanthus_passer-ONT-Illumina-Genome-Assembly/blob/master/1.%20Genome_Assembly_Ont_Illumina.md#trim-nanopore-adapters-and-filter-by-quality-and-length
#Filter to keep only longer but plausibly worse reads for assembly
gunzip -c ${porechop_out} | NanoFilt -q 3 -l 1000 | gzip > ${nanofilt_assemble}

#Filter to keep shorter but better reads for polishing
gunzip -c ${porechop_out} | NanoFilt -q 5 -l 500 | gzip > ${nanofilt_polish}

#Output number of reads etc. of the two assemblies
conda deactivate; module load anaconda3; source activate nanoplot
NanoPlot \
  -t ${SLURM_CPUS_PER_TASK} \
  --fastq ${nanofilt_assemble} \
  --plots hex dot \
  --title pc_nA \
  --outdir ${outDIR}/nanostat_pc_nA/

NanoPlot \
  -t ${SLURM_CPUS_PER_TASK} \
  --fastq ${nanofilt_polish} \
  --plots hex dot \
  --title pc_nP \
  --outdir ${outDIR}/nanostat_pc_nP/
