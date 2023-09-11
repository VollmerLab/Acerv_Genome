#!/bin/bash

#SBATCH --job-name=removeMTDNA
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#https://github.com/samtools/samtools/issues/1639
#genome_assembly=/work/vollmer/k2_nanopore_genome/k2_flyeHQ_polish.fa
#mitogenome_read_assembly=/scratch/j.selwyn/tmp/mitogenome/AC1_DV_GCCAAT-fp-fqscrn-repr/mt_assembly/AC1_DV_GCCAAT-fp-fqscrn-repr.megahit.result/AC1_DV_GCCAAT-fp-fqscrn-repr.megahit.mitogenome.fa

genome_assembly=${1}
mitogenome_read_assembly=${2}

outName=${genome_assembly%.*}_mtRemove.${genome_assembly##*.}

module load minimap2/2.17
module load samtools/1.10

cd $(mktemp -d)

#Find mitochondrial contigs
minimap2 -ax asm20 ${mitogenome_read_assembly} ${genome_assembly} -t ${SLURM_CPUS_PER_TASK} -2 > mapped.sam
#samtools view -b -f 4 ${outdir}/mapped.sam > genome.bam
samtools view -b -F 4 mapped.sam > mitochondrion.bam
mtContigs=$(samtools view mitochondrion.bam | cut -f1 | uniq)
printf "${mtContigs}\n" > mitoContigs.id

#Remove those contigs from genome
module load anaconda3; source activate seqtk

grep -o '^>.*' ${genome_assembly} | sed -e 's/>//g' > all_contigs.id
keep_ids=($(awk '{print $1}' all_contigs.id | grep -v -f mitoContigs.id))
printf "%s\n" "${keep_ids[@]}" > nuclear.id
seqtk subseq ${genome_assembly} nuclear.id > ${outName}
