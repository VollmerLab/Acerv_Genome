#!/bin/bash

#SBATCH --job-name=pruge_dups
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#outdir=/scratch/j.selwyn/tmp/purge_dupes
#in_genome=/scratch/j.selwyn/tmp/flye/k2HybridAssembly_longreads_flyeHQ.fasta
#longreads=/scratch/j.selwyn/tmp/raw_longreads/k2HybridAssembly_nanopore_input.fastq.gz

outdir=${1}
in_genome=${2}
longreads=${3}

module load minimap2/2.17
module load anaconda3; source activate purge_dups

genome_base=$(echo ${in_genome} | sed -e 's!.*/!!')

mkdir -p ${outdir}
ln -s ${in_genome} ${outdir}
cd ${outdir}
in_genome=$(pwd)/$(ls *fasta)

## Step 1: get coverage across genome
minimap2 -t ${SLURM_CPUS_PER_TASK} -xmap-ont ${in_genome} ${longreads} | gzip -c - > ${outdir}/genome.paf.gz
pbcstat -O ${outdir} ${outdir}/genome.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts ${outdir}/PB.stat > ${outdir}/cutoffs 2> ${outdir}/calcults.log

## Step 2: Split an assembly and do a self-self alignment
split_fa ${in_genome} > ${in_genome}.split
minimap2 -t ${SLURM_CPUS_PER_TASK} -xasm5 -DP ${in_genome}.split ${in_genome}.split | gzip -c - > ${in_genome}.split.self.paf.gz

## Step 3: Purge Haplotigs & Dupes
purge_dups -2 -T ${outdir}/cutoffs -c ${outdir}/PB.base.cov ${in_genome}.split.self.paf.gz > ${outdir}/dups.bed 2> ${outdir}/purge_dups.log

## Step 4: Remove Duplicates & Haplotigs
get_seqs -e ${outdir}/dups.bed ${in_genome}

mv ${outdir}/purged.fa ${outdir}/${genome_base%%.*}_purgeDup.${genome_base##*.}

# Spawn assembly check job
scriptDir=/work/vollmer/software/jds_scripts
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_base%%.*}_purgeDup.${genome_base##*.}
