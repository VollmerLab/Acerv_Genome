#!/bin/bash

#SBATCH --job-name=cafe
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
LSD2=/work/vollmer/software/lsd2.sif
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
species_tree=${2}
time_calibration=${3}
loci_dir=${4}
gene_families=${5}
max_diff=${6}

#outdir=/scratch/j.selwyn/phylogenetics/time_calibration_cafe
#species_tree=/scratch/j.selwyn/phylogenetics/orthofinder/Results_Mar29/Species_Tree/SpeciesTree_rooted.txt
#time_calibration=/scratch/j.selwyn/phylogenetics/fossil_estimates.txt
#loci_dir=/scratch/j.selwyn/phylogenetics/orthofinder/Results_Mar29/Orthogroup_Sequences
#gene_families=/scratch/j.selwyn/phylogenetics/orthofinder/Results_Mar29/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
#max_diff=100

mkdir -p ${outdir}
cd ${outdir}
ln -s ${species_tree} ${outdir}/tree.orthofinder

#Get concatenated length of all orthogroup sequences
seqLength=$(bash ${scriptDir}/runRscript.slurm \
  ${scriptDir}/r_utils/orthogroup_sequence_length.R \
    ${loci_dir})
    NUMBER=$()
seqLength=$(echo ${seqLength} | sed 's/[^0-9]*//g')

#read tree into R and then write.tree it back out to make proper newick format??
singularity exec --bind /work,/scratch,/tmp ${LSD2} \
  tree_convert \
    ${outdir}/tree.orthofinder

#need to figure out what -s means/is/hpw to figure out given orthofinder
singularity exec --bind /work,/scratch,/tmp ${LSD2} \
  lsd2 \
    -i ${outdir}/tree.orthofinder \
    -d ${time_calibration} \
    -z 0 \
    -s ${seqLength} \
    -f 1000

#Make tree ultrametric to deal with issues caused by rounding
singularity exec --bind /work,/scratch,/tmp ${LSD2} \
  tree_ultrametric \
    ${outdir}/tree.orthofinder.result.date.nexus

#Convert Orthofinder into proper format for cafe
bash ${scriptDir}/runRscript.slurm \
  ${scriptDir}/r_utils/cafe5_input_creation.R \
    ${gene_families} \
    ${outdir}/gene_families.txt \
    ${max_diff}

#Run Cafe analysis
module load anaconda3; source activate cafe
mkdir -p ${outdir}/cafeOut

cafe5 \
  --cores ${SLURM_CPUS_PER_TASK} \
  --output_prefix ${outdir}/cafeOut/errorModel \
  --infile ${outdir}/gene_families.txt \
  --tree ${outdir}/tree.orthofinder.result.date.nexus.ultrametric \
  --poisson \
  -e

# 2 - run lambda model with estimated error
cafe5 \
  --cores ${SLURM_CPUS_PER_TASK} \
  --output_prefix ${outdir}/cafeOut/errorModel \
  --infile ${outdir}/gene_families.txt \
  --tree ${outdir}/tree.orthofinder.result.date.nexus.ultrametric \
  --poisson \
  -e${outdir}/cafeOut/errorModel/Base_error_model.txt


for k in {1..10}; do
  # 1 - run gamma models with different numbers of bins - no error estimates though
  cafe5 \
    --cores ${SLURM_CPUS_PER_TASK} \
    -k ${k} \
    --output_prefix ${outdir}/cafeOut/K${k} \
    --infile ${outdir}/gene_families.txt \
    --tree ${outdir}/tree.orthofinder.result.date.nexus.ultrametric \
    --poisson
done
