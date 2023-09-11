#!/bin/bash

#SBATCH --job-name=blobFilter
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
BLOBTOOLS=/work/vollmer/software/blobtools2.sif
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
draft_genome=${2}
MINLENGTH=${3}
MINCOVERAGE=${4}

taxid=6130 #https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6130
protein_database=/work/vollmer/Databases/blobtools/uniprot/reference_proteomes.dmnd
taxdump_directory=/work/vollmer/Databases/blobtools/taxdump
BUSCO_database=/work/vollmer/Databases/busco_datasets/

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

MEMORY_PER_TASK=$(echo "scale=2 ; 0.8 * ${SLURM_MEM_PER_NODE} / ${SLURM_CPUS_PER_TASK}" | bc)
MEMORY_PER_TASK=${MEMORY_PER_TASK%.*}

mkdir -p ${outdir}
cd ${outdir}

## Run BUSCO to add to blob
module load anaconda3; source activate busco
busco \
  -i ${draft_genome} \
  -o busco \
  -m genome \
  --out_path=${outdir} \
  -c ${SLURM_CPUS_PER_TASK} \
  -l metazoa_odb10 \
  -f \
  --download_path=/work/vollmer/Databases/busco_datasets
conda deactivate

singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
  blobtools \
    add \
      --busco ${outdir}/busco/run_metazoa_odb10/full_table.tsv \
      --threads ${SLURM_CPUS_PER_TASK} \
      ${outdir}/${genome_name}_blob

## Interactive doesn't seem to work...
#singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
#   blobtools view --interactive ${outdir}/${genome_name}_blob

## Make plots of various metrics
#plot_options="blob cumulative snail"
#plot_outputs="svg png"
#for plotType in ${plot_options}
#do
#  for outType in ${plot_outputs}
#  do
#    export VIEWER=true; singularity run --bind /work,/scratch,/tmp ${BLOBTOOLS}
#    singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
#      blobtools \
#        view \
#          --host http://localhost:8080 \
#          --driver chromium \
#          --format png \
#          --out . \
#          --view cumulative \
#          ${outdir}/${genome_name}_blob
#  done
#done

## Get Summary Stats of Contigs
singularity exec --bind /work,/scratch,/tmp ${BLOBTOOLS} \
  blobtools \
    filter \
      --taxdump ${taxdump_directory} \
      --table prefilter_contig_info.tsv \
      --table-fields gc,length,merged_reads_cov,bestsumorder_phylum,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum_score \
      ${outdir}/${genome_name}_blob

## Make some Plots & Extract
module load R/4.0.3
Rscript ${scriptDir}/r_utils/blobtoolsFilter.R \
  ${outdir} \
  ${MINLENGTH} \
  ${MINCOVERAGE}

## Extract only contigs that passed filter (metazoan & length/coverage)
contigs_to_keep=$(cut -f2 ${outdir}/retained_contig_info.tsv)

source activate seqtk
touch ${outdir}/${genome_name}_blob.${draft_genome##*.}
for contig in ${contigs_to_keep}; do
  echo ${contig}
  printf "${contig}\n" | seqtk subseq ${draft_genome} - >> ${outdir}/${genome_name}_blob.${draft_genome##*.}
done

# Spawn assembly check job
conda deactivate
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_blob.${draft_genome##*.}
