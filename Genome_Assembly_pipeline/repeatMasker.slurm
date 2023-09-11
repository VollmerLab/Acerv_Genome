#!/bin/bash

#SBATCH --job-name=RepeatMasker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=150G

module load singularity
MAKER=/work/vollmer/software/maker.sif

#genome=/work/vollmer/k2_nanopore_genome/k2_flyeHQ_polish.fa
#outdir=/scratch/j.selwyn/tmp/repeatMaskerFull
#repeatModel=/scratch/j.selwyn/tmp/repeatModelerFull/k2_flyeHQ_polish-families.fa
#minLength_genome=1000000


outdir=${1}
genome=${2}
repeatModel=${3}

mkdir -p ${outdir}
cd ${outdir}

genome_name=$(echo ${genome##*/} | sed -e "s/\.f.*a$//")

#module load anaconda3; source activate bioawk
#bioawk -v LENGTH="${minLength_genome}" -c fastx '{ if(length($seq) > LENGTH) { print ">"$name; print $seq }}' \
#  ${genome} > ${outdir}/${genome_name}_${minLength_genome}.fa
#conda deactivate

#Repeat Masker
export BLAST_USAGE_REPORT=FALSE
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  RepeatMasker \
  -lib ${repeatModel} \
  -pa $(( SLURM_CPUS_PER_TASK / 4 )) \
  -engine ncbi \
  -dir ${outdir} \
  -s \
  ${outdir}/${genome_name}.fa


#Convert merged repeat masker to gff3 format
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  rmOutToGFF3.pl \
    ${outdir}/${genome_name}.fa.out > ${outdir}/${genome_name}.out.gff3

# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" ${outdir}/${genome_name}.out.gff3 \
  > ${outdir}/${genome_name}.complex.gff3
# reformat to work with MAKER
cat ${outdir}/${genome_name}.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
  > ${outdir}/${genome_name}.complex.reformat.gff3
