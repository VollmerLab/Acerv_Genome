#!/bin/bash

#SBATCH --job-name=allBUSCO
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

scriptDir=/work/vollmer/software/jds_scripts
module load anaconda3; source activate busco

outdir=${1}
proteinDir=${2}

#outdir=/scratch/j.selwyn/phylogeny_proteome3/buscos
#proteinDir=/scratch/j.selwyn/phylogeny_proteome3/proteins

mkdir -p ${outdir}
cd ${outdir}

all_species=$(ls ${proteinDir}/*fasta)

for species_proteins in ${all_species}; do
  busco_extension=$(echo ${species_proteins##*/} | sed "s/.fasta//g")

  echo ${species_proteins}
  echo ${busco_extension}

  busco \
    -i ${species_proteins} \
    -o busco \
    -m proteins \
    --out_path=${outdir}/${busco_extension} \
    -c ${SLURM_CPUS_PER_TASK} \
    -l metazoa_odb10 \
    -f \
    --download_path=/work/vollmer/Databases/busco_datasets \
    --update-data
done

conda deactivate
bash ${scriptDir}/runRscript.slurm \
  ${scriptDir}/r_utils/summarize_busco_dir.R \
    ${outdir}
