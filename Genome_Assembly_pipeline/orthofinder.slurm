#!/bin/bash

#SBATCH --job-name=OrthoFinder
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G

#Cluster to longest isoform!

module load anaconda3; source activate orthofinder

outdir=${1}
indir=${2}
tree_type=${3:-default}

#outdir=/scratch/j.selwyn/phylogeny/orthofinder
#indir=/scratch/j.selwyn/phylogeny/proteins
#tree_type=default


echo ${tree_type}

if [ ${tree_type} = msa ]; then
  echo "Running Orthofinder with IQTree will spawn two subsequent jobs"
  DAY=$(date -d "$D" '+%d')
  MONTH=$(date -d "$D" '+%b')
  workingDir=Results_${MONTH}${DAY}

  mkdir -p ${outdir}
  cd ${outdir}

  orthofinder \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${SLURM_CPUS_PER_TASK} \
    -fg ${indir} \
    -M msa

else
  echo "Running Orthofinder with Default Tree Builder"
  DAY=$(date -d "$D" '+%d')
  MONTH=$(date -d "$D" '+%b')
  workingDir=Results_${MONTH}${DAY}

  orthofinder \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${SLURM_CPUS_PER_TASK} \
    -o ${outdir} \
    -f ${indir}

  #After running orthofinder then spawn job to remake tree using MSA alignment to get bootstrap support
  #sbatch \
  #  --output=${outdir}/../slurm/OrthoFinder_iqtree_%j.output \
  #  /work/vollmer/software/jds_scripts/orthofinder.slurm \
  #    ${outdir}/../ortho_msa \
  #    ${outdir}/${workingDir} \
  #    ${timeCalibration} \
  #    msa

  #sbatch \
  #  --output=${outdir}/../slurm/time_calibration_cafe_%j.output \
  #  /work/vollmer/software/jds_scripts/runCAFE.slurm \
  #    ${outdir}/../time_calibration_cafe \
  #    ${outdir}/${workingDir}/Species_Tree/SpeciesTree_rooted.txt \
  #    ${timeCalibration} \
  #    ${outdir}/${workingDir}/Orthogroup_Sequences \
  #    ${outdir}/${workingDir}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv \
  #    100
fi
