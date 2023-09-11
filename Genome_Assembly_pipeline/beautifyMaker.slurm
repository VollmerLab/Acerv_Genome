#!/bin/bash

#SBATCH --job-name=rename_maker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#https://www.bioinformatics.uni-muenster.de/publication_data/P.californicus_annotation/functional_prediction.hbi

outdir=${1}
maker_gff=${2}
maker_transcripts=${3}
maker_proteins=${4}
genome=${5}

module load singularity
MAKER=/work/vollmer/software/maker.sif
scriptDir=/work/vollmer/software/jds_scripts

cd ${outdir}
base_name=$(echo ${maker_gff##*/} | sed -e "s/\.gff$//")

#Rename everything to how NCBI suggests
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  maker_map_ids \
    --prefix Acer_ \
    --justify 8 \
    ${maker_gff} > ${outdir}/${base_name}.id.map

#Copy original files to be renamed
cp ${maker_gff} ${outdir}/.
cp ${maker_transcripts} ${outdir}/.
cp ${maker_proteins} ${outdir}/.

#Rename with new names
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  map_gff_ids \
    ${outdir}/${base_name}.id.map \
    ${outdir}/${base_name}.gff

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  map_fasta_ids \
    ${outdir}/${base_name}.id.map \
    ${outdir}/${base_name}.all.maker.transcripts.fasta

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  map_fasta_ids \
    ${outdir}/${base_name}.id.map \
    ${outdir}/${base_name}.all.maker.proteins.fasta

#Remove rows in gff containing exclusively "###"
sed -i '/^###$/d' ${outdir}/${base_name}.gff

#Spawn R Job to produce cleaner easier to use structural annotation files and decoder file between Trinity transcripts and Maker genes
sbatch \
  --cpus-per-task=${SLURM_CPUS_PER_TASK} \
  --job-name=postProcessStructural \
  --output=${outdir}/slurm/postProcessStructural_%j.output \
  ${scriptDir}/runRscript.slurm \
    ${scriptDir}/r_utils/postprocess_structural_annotations.R \
      ${outdir}/${base_name}.gff \
      ${genome} \
      ${SLURM_CPUS_PER_TASK}
