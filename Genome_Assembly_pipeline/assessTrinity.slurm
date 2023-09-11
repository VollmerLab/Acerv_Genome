#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=trinityStats
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
#https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-without-a-genome.html#gsc.tab=0

#transcriptome_assembly=/scratch/j.selwyn/genome_annotation/flyeHQ_16-July-2022/trinity/k2_flyeHQ_polish_trinity-gg.fasta
#rawSequence=/scratch/j.selwyn/genome_annotation/flyeHQ_16-July-2022/repair/*gz

transcriptome_assembly=${1}
#; shift
#rawSequence="$@"
#echo ${rawSequence}
#samples_file=${2} #see below for example.

#      --samples_file <string>    tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.


working_dir=${transcriptome_assembly%/*}
genome_name=$(echo ${transcriptome_assembly##*/} | sed -e "s/_trinity-gg.fasta$//")

cd ${working_dir}

#cp ${samples_file} ${working_dir}/.

module load samtools
module load singularity

TRINITY=/work/vollmer/software/trinity.sif
GMAP=/work/vollmer/software/gmap.sif


#Trinity Contiguity Stats
singularity exec --bind /work,/scratch,/tmp ${TRINITY} \
   /usr/local/bin/util/TrinityStats.pl \
   ${transcriptome_assembly} > ${working_dir}/Trinity_assembly.metrics

## BUSCO
module load anaconda3; source activate busco

busco \
  -i ${transcriptome_assembly} \
  -o busco \
  -m transcriptome \
  --out_path=${working_dir}/ \
  -c ${SLURM_CPUS_PER_TASK} \
  -l metazoa_odb10 \
  -f \
  --download_path=/work/vollmer/Databases/busco_datasets \
  --update-data
conda deactivate
