#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=indexGenome
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

# Script that takes in a genome and indexes it as either DNA (with bwa index) or RNA (with gmap_build)

outdir=${1}
ref_genome=${2}
nucleotide=${3}
ref_gff3=${4} #must include when using 'RNA_hisat'


if [ ${nucleotide} != "DNA" ] && [ ${nucleotide} != "RNA" ] && [ ${nucleotide} != "RNA_hisat" ]
then
  echo "Nucleotide (\$3) must be one of DNA or RNA or RNA_hisat"
  exit 0
else
  echo "Running with Nucleotide ${nucleotide}"
fi


mkdir -p ${outdir}
cd ${outdir} #move to outdir for if temp files are made
ln -s ${ref_genome} ${outdir}/reference.fa

module load bwa; module load samtools/1.10

if [ ${nucleotide} == "DNA" ]
then
  bwa index ${outdir}/reference.fa

elif [ ${nucleotide} == "RNA" ]
then
  #Build Genome database
  module load singularity
  CONTAINER=/work/vollmer/software/gmap.sif

  singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
    gmap_build \
      --dir=${outdir}/genomeDB \
      -d=genomeDB \
      -t ${SLURM_CPUS_PER_TASK} \
      ${outdir}/reference.fa

else
  module load singularity
  AGAT=/work/vollmer/software/agat.sif
  ref_gtf=$(echo ${ref_gff3} | sed -E 's/gff(3)?/gtf/g')

  singularity exec --bind /work,/scratch,/tmp ${AGAT} \
    agat_convert_sp_gff2gtf.pl \
    --gff ${ref_gff3} \
    -o ${ref_gtf}

  module load anaconda3; source activate hisat2
  hisat2_extract_splice_sites.py ${ref_gtf} > ${outdir}/splicesites.tsv
  hisat2_extract_exons.py -v ${ref_gtf} > ${outdir}/exons.tsv

  hisat2-build \
    -p ${SLURM_CPUS_PER_TASK} \
    --ss ${outdir}/splicesites.tsv \
    --exon ${outdir}/exons.tsv \
    ${outdir}/reference.fa \
    ${outdir}/reference_index
fi

grep -o '^>.*' ${outdir}/reference.fa | sed -e 's/>//g' > ${outdir}/reference.contigs
samtools faidx ${outdir}/reference.fa
