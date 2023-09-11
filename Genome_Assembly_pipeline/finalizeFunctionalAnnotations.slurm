#!/bin/bash

#SBATCH --job-name=FinalInterPro
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

#base_name=r3

module load singularity
# module load R/4.0.3
MAKER=/work/vollmer/software/maker.sif
RSTUDIO_IMAGE="/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif"
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
genome=${2}
base_name=${3}

cd ${outdir}

#Gather BLAST/uniprot
cat ${outdir}/blast_array/*/*blastp > ${outdir}/${base_name}.blast.out
#cp ${outdir}/${base_name}.blast.out ${outdir}/${base_name}.blast_rename.out
#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  map_data_ids \
#    ${outdir}/${base_name}.id.map \
#    ${outdir}/${base_name}.blast_rename.out

#Gather InterProt
cat ${outdir}/interprot_array/*/protein_chunks.fasta.tsv > ${outdir}/${base_name}.interprot.tsv
#cp ${outdir}/${base_name}.interprot.tsv ${outdir}/${base_name}.interprot_rename.tsv
#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  map_data_ids \
#  ${outdir}/${base_name}.id.map \
#  ${outdir}/${base_name}.interprot_rename.tsv

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  gff3_merge \
  -o ${outdir}/${base_name}.interprot.gff3 \
  ${outdir}/interprot_array/*/protein_chunks.fasta.gff3

#Gather KOFAM
cat ${outdir}/kofamscan_array/*/*kofamscan > ${outdir}/${base_name}.kofamscan.out

#Gather EnTAP
output_types=$(ls ${outdir}/entap_array/0/final_results | sed -e 's/.*\(\..*\..*$\)/\1/g' | sort -u)
for output in ${output_types}; do
  echo ${output}
  cat ${outdir}/entap_array/*/final_results/${output} > ${outdir}/TMP${output}

  #remove rows that match first row excluding first row - from tsv files not fasta style
  if [[ ${output} == *"tsv" ]]; then
    header_line=$(cat ${outdir}/TMP${output} | head -n 1)
    cat ${outdir}/TMP${output} | sed "s/${header_line}//g" | sed '/^$/d' | sed "1 i ${header_line}" | sed 's/\t$//g' >\
      ${outdir}/${base_name}.${output}
    rm ${outdir}/TMP${output}
  else
    mv ${outdir}/TMP${output} ${outdir}/${base_name}.${output}
  fi
  rename final entap ${outdir}/${base_name}.${output}
done


#Clean up folders
# rm -rf ${outdir}/blast_array
# rm -rf ${outdir}/interprot_array
# rm -rf ${outdir}/kofamscan_array

#R script to merge into rds output
singularity run -B "/home:/home,/scratch:/scratch,/work:/work" ${RSTUDIO_IMAGE} \
  Rscript \
    ${scriptDir}/r_utils/merge_functional_annotations.R \
      ${outdir} \
      ${base_name}.blast.out \
      ${base_name}.interprot.tsv \
      ${base_name}.kofamscan.out \
      ${base_name}.entap_annotations_lvl0.tsv

#Map functional from BLAST to gff & fastas
#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  maker_functional_gff \
#    ${genome} \
#    ${outdir}/${base_name}.blast_rename.out \
#    ${outdir}/${base_name}.gff > ${outdir}/${base_name}.putative_function.gff

#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  maker_functional_fasta \
#    ${genome} \
#    ${outdir}/${base_name}.blast_rename.out \
#    ${outdir}/${base_name}.all.maker.proteins.fasta > ${outdir}/${base_name}.all.maker.proteins.putative_function.fasta

#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  maker_functional_fasta \
#    ${genome} \
#    ${outdir}/${base_name}.blast_rename.out \
#    ${outdir}/${base_name}.all.maker.transcripts.fasta > ${outdir}/${base_name}.all.maker.transcripts.putative_function.fasta

#Add interprotscan outputs
#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  ipr_update_gff \
#    ${outdir}/${base_name}.putative_function.gff \
#    ${outdir}/${base_name}.interprot_rename.tsv > ${outdir}/${base_name}.putative_domains.gff

#singularity exec --bind /work,/scratch,/tmp ${MAKER} \
#  iprscan2gff3 \
#    ${outdir}/${base_name}.interprot_rename.tsv \
#    ${outdir}/${base_name}.gff > ${outdir}/${base_name}_visible_iprscan_domains.gff
