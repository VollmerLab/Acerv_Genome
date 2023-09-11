#!/bin/bash

#SBATCH --job-name=finalizeMaker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

module load singularity
MAKER=/work/vollmer/software/maker.sif

outdir=${1}
annotation_round=${2}

cat ${outdir}/maker_array/datastore_index_*.log > ${outdir}/maker_array/r${annotation_round}_datastore_index.log


singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  gff3_merge -s -d ${outdir}/maker_array/r${annotation_round}_datastore_index.log > \
  ${outdir}/r${annotation_round}.gff

#not sure what this should be doing...
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  fasta_merge \
    -d ${outdir}/maker_array/r${annotation_round}_datastore_index.log \
    -o ${outdir}/r${annotation_round}

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  gff3_merge -n -s -d ${outdir}/maker_array/r${annotation_round}_datastore_index.log > \
  ${outdir}/r${annotation_round}.noseq.gff

#rm -rf ${outdir}/maker_array

#Filter out poor AED fits to be used in next MAKER rounds
#module load anaconda3; source activate agat
#agat_sp_select_feature_by_attribute_value.pl \
#  --gff ${outdir}/r${annotation_round}.gff \
#  --value 1 \
#  -t "=" [ --output outfile ]

#
#singularity exec --bind /work,/scratch,/tmp ${MAKER} quality_filter.pl -a 0.5 r3.gff > r3_filter.gff
