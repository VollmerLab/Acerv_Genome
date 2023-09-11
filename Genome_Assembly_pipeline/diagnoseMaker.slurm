#!/bin/bash

#SBATCH --job-name=diagnoseMaker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
module load R/4.0.3
MAKER=/work/vollmer/software/maker.sif
scriptDir=/work/vollmer/software/jds_scripts

#full_gff=/scratch/j.selwyn/tmp/maker_r1/r1.gff
full_gff=${1}
full_transcripts=${2}

outdir=${full_gff%/*}
cd ${outdir}

cat ${full_gff} | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' > ${outdir}/nGene_lengthGene.txt

singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  AED_cdf_generator.pl -b 0.025 ${full_gff} > ${outdir}/aed_cdf.dat

Rscript ${scriptDir}/r_utils/aed_plot.R ${outdir}/aed_cdf.dat

#
module load anaconda3; source activate busco
export AUGUSTUS_CONFIG_PATH=${outdir}/../augustus_config/config
busco \
  -i ${full_transcripts} \
  -o annotation_eval \
  -m transcriptome \
  --out_path=${outdir}/ \
  -c ${SLURM_CPUS_PER_TASK} \
  --augustus_species acropora_cervicornis \
  --augustus_parameters='--progress=true' \
  -l metazoa_odb10 \
  -f \
  --download_path=/work/vollmer/Databases/busco_datasets \
  --update-data
