#!/bin/bash

#SBATCH --job-name=Trinity
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G

mapped_reference=${1}
maxIntron=${2}

working_dir=${mapped_reference%/*}
genome_name=$(echo ${mapped_reference##*/} | sed -e "s/_sorted.bam$//")

TOTAL_MEM=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE}" | bc)
TOTAL_MEM_GB=$(echo "scale=0 ; ${TOTAL_MEM}/100" | bc)

MEMORY_PER_TASK=$(echo "scale=2 ; ${TOTAL_MEM} / ${SLURM_CPUS_PER_TASK}" | bc)

cd ${working_dir}

module load samtools
module load singularity

TRINITY=/work/vollmer/software/trinity.sif
GMAP=/work/vollmer/software/gmap.sif

#Start Trinity Run
singularity exec --bind /work,/scratch,/tmp ${TRINITY} \
  Trinity \
    --CPU ${SLURM_CPUS_PER_TASK} \
    --max_memory ${TOTAL_MEM_GB}G \
    --genome_guided_bam ${mapped_reference} \
    --genome_guided_max_intron ${maxIntron} \
    --output ${working_dir}/trinity_out

cp ${working_dir}/trinity_out/Trinity-GG.fasta ${working_dir}/${genome_name}_trinity-gg.fasta

singularity exec --bind /work,/scratch,/tmp ${TRINITY} \
  /usr/local/bin/util/misc/cdna_fasta_file_to_transcript_gtf.pl \
    ${working_dir}/${genome_name}_trinity-gg.fasta > ${working_dir}/${genome_name}_trinity-gg.gtf

module load anaconda3; source activate gffread
gffread ${working_dir}/${genome_name}_trinity-gg.gtf > ${working_dir}/${genome_name}_trinity-gg.gff3
conda deactivate

#Map Trinity Transcripts to Genome - i probably accidentally deleted this once while debugging code
if [ ! -d ${working_dir}/genomeDB ]
then
  #Build Genome database
  singularity exec --bind /work,/scratch,/tmp ${GMAP} \
    gmap_build \
      --dir=${working_dir}/genomeDB \
      -d=genomeDB \
      -t ${SLURM_CPUS_PER_TASK} \
      ${working_dir}/reference.fa
fi

singularity exec --bind /work,/scratch,/tmp ${GMAP} \
  gmap \
    -D ${working_dir}/genomeDB \
    -d genomeDB \
    -f samse \
    -n 0 \
    -x 50 \
    -t ${SLURM_CPUS_PER_TASK} \
    -B 4 \
    ${working_dir}/${genome_name}_trinity-gg.fasta > ${working_dir}/${genome_name}_trinity-gg.sam

samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${working_dir}/${genome_name}_trinity-gg.sam > ${working_dir}/transcriptome_to_genome_mapping_stats
