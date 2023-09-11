#!/bin/bash

#SBATCH --job-name=finalizeSC
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G

module load singularity
module load anaconda3; source activate mummer
CONTAINER=/work/vollmer/software/finishingtool.sif
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
draft_genome=${2}
long_polishing_reads=${3}

# ; shift 3
# short_polishing_reads=( "$@" )
# short_read_only_assembly=${4} - use short read assembly for gap filling too??

genome_name=$(echo ${draft_genome##*/} | sed -e "s/\.f.*a//")

mkdir -p ${outdir}
cd ${outdir}


if [ -f "${outdir}/improved2.fasta" ]; then
    echo "Intermediate File exists."
    echo "Picking up from: improved2.fasta"

    singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
      python \
        /opt/finishingTool/finisherSC.py \
          -par ${SLURM_CPUS_PER_TASK} \
          -f True \
          -p ${outdir}/improved2.fasta \
          ${outdir} \
          /work/vollmer/software/.conda/envs/mummer/bin
else
  echo "Intermediate File does not exist"
  echo "Initiating FinisherSC"
  ln -s ${draft_genome} ${outdir}/contigs.fasta
  ln -s ${long_polishing_reads} ${outdir}/raw_reads.fasta

  singularity exec --bind /work,/scratch,/tmp ${CONTAINER} \
    python \
      /opt/finishingTool/finisherSC.py \
        -par ${SLURM_CPUS_PER_TASK} \
        -f True \
        ${outdir} \
        /work/vollmer/software/.conda/envs/mummer/bin
fi

mv ${outdir}/improved3.fasta ${outdir}/${genome_name}_finishSC.${draft_genome##*.}

# Spawn assembly check job
conda deactivate
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_finishSC.${draft_genome##*.}
