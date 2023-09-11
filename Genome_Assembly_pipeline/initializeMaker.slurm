#!/bin/bash

#SBATCH --job-name=initMaker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

#genome=/work/vollmer/k2_nanopore_genome/k2_flyeHQ_polish.fa
#transcriptome=/scratch/j.selwyn/genome_annotation/flyeHQ_18-July-2022/trinity/k2_flyeHQ_polish_trinity-gg.fasta
#outdir=/scratch/j.selwyn/tmp2
#minLength_genome=1000000
#minLength_transcriptome=10000

outdir=${1}
genome=${2}
transcriptome=${3}
minLength_genome=${4}
minLength_transcriptome=${5}

mkdir -p ${outdir}
cd ${outdir}

genome_name=$(echo ${genome##*/} | sed -e "s/\.f.*a$//")
transcriptome_name=$(echo ${transcriptome##*/} | sed -e "s/\.f.*a$//")

#Shrink genome for testing purposes - may leave in too
module load anaconda3; source activate bioawk
bioawk -v LENGTH="${minLength_genome}" -c fastx '{ if(length($seq) > LENGTH) { print ">"$name; print $seq }}' \
  ${genome} > ${outdir}/${genome_name}_${minLength_genome}.fa

module load samtools
samtools faidx ${outdir}/${genome_name}_${minLength_genome}.fa

module load anaconda3; source activate assembly_stats
echo "Original Reference Genome Stats"
assembly-stats ${genome}
echo
echo
echo
echo "Filtered Reference Genome Stats"
assembly-stats ${outdir}/${genome_name}_${minLength_genome}.fa
conda deactivate

module load anaconda3; source activate bioawk
bioawk -v LENGTH="${minLength_transcriptome}" -c fastx '{ if(length($seq) > LENGTH) { print ">"$name; print $seq }}' \
  ${transcriptome} > ${outdir}/${transcriptome_name}_${minLength_transcriptome}.fa
conda deactivate

#Make lists of all contigs in genome and transcriptome
grep -o '^>.*' ${outdir}/${genome_name}_${minLength_genome}.fa | sed -e 's/>//g' > ${outdir}/${genome_name}_${minLength_genome}.contig
readarray -t genome_contigs < ${outdir}/${genome_name}_${minLength_genome}.contig

grep -o '^>.*' ${outdir}/${transcriptome_name}_${minLength_transcriptome}.fa | sed -e 's/>//g' > ${outdir}/${transcriptome_name}_${minLength_transcriptome}.contig
readarray -t transcriptome_contigs < ${outdir}/${transcriptome_name}_${minLength_transcriptome}.contig

n_genome=${#genome_contigs[@]}
n_transcriptome=${#transcriptome_contigs[@]}

#Make all combinations of contigs
#touch ${outdir}/genome.contig
#for gen_cont in ${genome_contigs[@]}; do
  #echo "Starting ${gen_cont}"
#  yes ${gen_cont} | head -n ${n_transcriptome} >> ${outdir}/genome.contig
#done

#touch ${outdir}/transcriptome.contig
#for (( i=1; i<=${n_genome}; i++ )); do
	#echo -n "$i "
#  printf "%s\n" "${transcriptome_contigs[@]}"  >> ${outdir}/transcriptome.contig
#done
