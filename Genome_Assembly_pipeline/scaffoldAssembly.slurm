#!/bin/bash

#SBATCH --job-name=scaffold
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load anaconda3; source activate longstitch
scriptDir=/work/vollmer/software/jds_scripts

outdir=${1}
draft_genome=${2}
polishing_reads=${3}

algo_options="tigmint-ntLink-arks"
K_choices="24 32 40"
W_choices="100 250 500"

tmp=${draft_genome##*/}
genome_name=${tmp%%\.*}

mkdir -p ${outdir}
#cd $(mktemp -d)
cd ${outdir}

#ln -s ${draft_genome} ./draft-assembly.fa
#ln -s ${polishing_reads} ./reads.fq.gz

## Grid Search to find best K/W settings
#Try just ntLink-arks i.e. no breaking contigs - currently no error corection/breaking of Misassembly
z=0; pids=()
#autospan=~10 (i think) #decreasing leads to bigger N50 and fewer bonus configs (i.e. fewer breaks and/or more repairs of breaks)
#dist=~2270 (i think) # increasing leads to bigger N50 and fewer bonus configs (i.e. fewer breaks and/or more repairs of breaks)
for K in ${K_choices}
do
  for W in ${W_choices}
  do
    mkdir -p ${outdir}/array_search/${z}
    cd ${outdir}/array_search/${z}

    ln -s ${draft_genome} ./draft-assembly.fa
    ln -s ${polishing_reads} ./reads.fq.gz

    if [ ${z} = "0" ]
    then
      longstitch \
        ${algo_options} \
          draft=draft-assembly \
          reads=reads \
          G=3.5e8 \
          z=100 \
          span=2 \
          dist=50000 \
          t=${SLURM_CPUS_PER_TASK} \
          gap_fill=True \
          rounds=5 \
          k_ntLink=${K} \
          w=${W}
    else
      longstitch \
        ${algo_options} \
          draft=draft-assembly \
          reads=reads \
          G=3.5e8 \
          z=100 \
          span=2 \
          dist=50000 \
          t=3 \
          gap_fill=True \
          rounds=5 \
          k_ntLink=${K} \
          w=${W} &
    fi
    z=$((z+1))
    pids+=( "$!" )
    cd ${outdir}
  done
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

## Calc N50 for each assembly
conda deactivate; module load anaconda3; source activate assembly_stats
z=0; all_scaffolds=(); n50=(); bestN50=1
for K in ${K_choices}
do
  for W in ${W_choices}
  do
    echo "K = ${K}; W = ${W}"
    all_scaffolds+=( "${outdir}/array_search/${z}/draft-assembly.k${K}.w${W}.${algo_options}.longstitch-scaffolds.fa" )

    stat_output=$(assembly-stats ${all_scaffolds[${z}]})

    tmpn50=$(echo ${stat_output} | grep -o -E "N50 = [0-9]+" | grep -o -E "[0-9]+$")
    n50+=( "${tmpn50}" )
    echo "N50 = ${n50[${z}]}"
    echo ""
    echo ""
    z=$((z+1))
  done
done

## Choose assembly with largest N50
for i in "${!n50[@]}"; do
  echo ${i}
        if ((n50[i]>arr1_max)); then
                arr1_max=${n50[i]}
                max_indexes=($i)
        elif [[ "${n50[i]}" == "$arr1_max" ]]; then
                max_indexes+=($i)
        fi
done

#mv $(readlink -f ${all_scaffolds[${max_indexes}]}) ${outdir}/${genome_name}_longstitch.${draft_genome##*.}

#Rename Contigs so Pilon doesn't get mad
conda deactivate; module load anaconda3; source activate bbmap
rename.sh in=$(readlink -f ${all_scaffolds[${max_indexes}]}) out=${outdir}/${genome_name}_longstitch.${draft_genome##*.} prefix=Acerv_scaffold

# Spawn assembly check job
conda deactivate
sbatch \
  --output=${outdir}/assembly_assessment_%j.output \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}_longstitch.${draft_genome##*.}
