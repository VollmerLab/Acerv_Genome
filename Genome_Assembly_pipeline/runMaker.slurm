#!/bin/bash

#SBATCH --job-name=Maker
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

module load singularity
MAKER=/work/vollmer/software/maker.sif
ARRAY_ID=${SLURM_ARRAY_TASK_ID}
#genome=/scratch/j.selwyn/tmp2/k2_flyeHQ_polish_1000000.fa
#transcriptome=/scratch/j.selwyn/tmp2/k2_flyeHQ_polish_trinity-gg_10000.fa
#outdir=/scratch/j.selwyn/tmp2

outdir=${1}
genome=${2}
transcriptome=${3}
protein_database=${4}
array_size=${5}
maker_round=${6}
seq_gff=${7} #if round 1 then unneeded
snap_train=${8} #if round 1 then unneeded
gm_train=${9} #if round 1 then unneeded

mkdir -p ${outdir}
cd ${outdir}

genome_name=$(echo ${genome##*/} | sed -e "s/\.f.*a$//")
transcriptome_name=$(echo ${transcriptome##*/} | sed -e "s/\.f.*a$//")

repeat_mask_complex=${outdir}/../${genome_name}.complex.reformat.gff3

split_jobs () {
  #https://www.geeksforgeeks.org/split-the-number-into-n-parts-such-that-difference-between-the-smallest-and-the-largest-part-is-minimum/
  X=$1
  N=$2

  out=""

  if [ $((X % N)) -eq 0 ]
    then
      for ((i=1; i<=N; i++))
      do
        out="${out} $((X / N))"
      done

    else
      zp=$((N - (X % N)))
      pp=$((X / N))

      for ((i=1; i<=N; i++))
      do
        if [ ${i} -gt ${zp} ]
        then
          out="${out} $((pp + 1))"
        else
          out="${out} ${pp}"
        fi
      done
    fi
  echo "${out}"
}

#readarray -t genome_contigs < ${outdir}/genome.contig
#readarray -t transcriptome_contigs < ${outdir}/transcriptome.contig

readarray -t genome_contigs < ${outdir}/../${genome_name}.contig
#something here to pick min of array size or number of contigs total
n_contig=${#genome_contigs[@]}

array_size=$(( n_contig < array_size ? n_contig : array_size ))

contig_bins=( $(split_jobs ${#genome_contigs[@]} ${array_size}) )

starting_index=$(echo "0 ${contig_bins[@]:0:$((ARRAY_ID))} 0" | sed -e 's/  */+/g' | bc)

genome_contig_names=( ${genome_contigs[@]:${starting_index}:${contig_bins[${ARRAY_ID}]}} )
#transcriptome_contig_names=( ${transcriptome_contigs[@]:${starting_index}:${contig_bins[${ARRAY_ID}]}} )

echo "Running maker on contigs: ${genome_contig_names[@]}"

tmpDir=${outdir}/maker_array/${ARRAY_ID}
mkdir -p ${tmpDir}
cd ${tmpDir}

#Subset the genome to the relevant portions
module load anaconda3; source activate seqtk
#printf '%s\n' "${genome_contig_names[@]}" | seqtk subseq ${genome} - > ${tmpDir}/genome_contigs.fasta
touch ${tmpDir}/genome_contigs.fasta
for contig in ${genome_contig_names[@]}; do
  echo ${contig}
  printf "${contig}\n" | seqtk subseq ${genome} - >> ${tmpDir}/genome_contigs.fasta
done
conda deactivate

#Make initial config files and adjust to specific genome & transcriptome chunks
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  maker \
  -CTL

#sed -i "s|^genome\=|genome\=${tmpDir}/genome_contigs.fasta|g" maker_opts.ctl
sed -i "s|^est\=|est\=${transcriptome}|g" maker_opts.ctl
sed -i "s|^model_org\=all|model_org\=simple|g" maker_opts.ctl
sed -i "s|^est2genome\=0|est2genome\=1|g" maker_opts.ctl
sed -i "s|^clean_up\=0|clean_up\=1|g" maker_opts.ctl
sed -i "s|^max_dna_len\=100000|max_dna_len\=10000000|g" maker_opts.ctl
sed -i "s|^rm_gff\=|rm_gff\=${repeat_mask_complex}|g" maker_opts.ctl
sed -i "s|^protein2genome\=0|protein2genome\=1|g" maker_opts.ctl
sed -i "s|^protein\=|protein\=${protein_database}|g" maker_opts.ctl
sed -i "s|^run_evm\=0|run_evm\=1|g" maker_opts.ctl
sed -i "s|^pred_stats\=0|pred_stats\=1|g" maker_opts.ctl
#sed -i "s|^min_contig\=1|min_contig\=10000|g" maker_opts.ctl


if [[ ${maker_round} -gt 1 ]]
then
  export AUGUSTUS_CONFIG_PATH=${outdir}/../augustus_config/config

  # transcript alignments
  #awk '{ if ($2 == "est2genome") print $0 }' ${noseq_gff} > ${tmpDir}/all.maker.est2genome.gff
  # protein alignments
  #awk '{ if ($2 == "protein2genome") print $0 }' ${noseq_gff} > ${tmpDir}/all.maker.protein2genome.gff
  # repeat alignments
  #awk '{ if ($2 ~ "repeat") print $0 }' ${noseq_gff} > ${tmpDir}/all.maker.repeats.gff

  sed -i "s|^model_org\=all|model_org\=|g" maker_opts.ctl
  sed -i "s|^est\=${transcriptome}|est\=|g" maker_opts.ctl
  sed -i "s|^protein\=${protein_database}|protein\=|g" maker_opts.ctl
  sed -i "s|^repeat_protein\=/usr/local/maker/data/te_proteins.fasta|repeat_protein\=|g" maker_opts.ctl
  sed -i "s|^rm_gff\=${repeat_mask_complex}|rm_gff\=|g" maker_opts.ctl

  sed -i "s|^maker_gff\=|maker_gff\=${seq_gff}|g" maker_opts.ctl
  sed -i "s|^est_pass\=0|est_pass\=1|g" maker_opts.ctl
  sed -i "s|^protein_pass\=0|protein_pass\=1|g" maker_opts.ctl
  sed -i "s|^rm_pass\=0|rm_pass\=1|g" maker_opts.ctl

  sed -i "s|^snaphmm\=|snaphmm\=${snap_train}|g" maker_opts.ctl
  sed -i "s|^augustus_species\=|augustus_species\=acropora_cervicornis|g" maker_opts.ctl
  sed -i "s|^gmhmm\=|gmhmm\=${gm_train}|g" maker_opts.ctl
  sed -i "s|^est2genome\=1|est2genome\=0|g" maker_opts.ctl
  sed -i "s|^protein2genome\=1|protein2genome\=0|g" maker_opts.ctl
  sed -i "s|^trna\=0|trna\=1|g" maker_opts.ctl

  #https://github.com/sujaikumar/assemblage/blob/master/README-annotation.md
  #sed -i "s|^min_protein\=0|min_protein\=30|g" maker_opts.ctl
  #sed -i "s|^alt_splice\=0|alt_splice\=1|g" maker_opts.ctl
  #sed -i "s|^keep_preds\=0|keep_preds\=1|g" maker_opts.ctl #https://groups.google.com/g/maker-devel/c/VaoXWlGHOjs
  #sed -i "s|^split_hit\=10000|split_hit\=4000|g" maker_opts.ctl
  #sed -i "s|^single_exon\=0|single_exon\=1|g" maker_opts.ctl
  #sed -i "s|^single_length\=250|single_length\=250|g" maker_opts.ctl
fi

#if [[ ${maker_round} -gt 2 ]]
#then
  #sed -i "s|^pred_pass\=0|pred_pass\=1|g" maker_opts.ctl
  #sed -i "s|^model_pass\=0|model_pass\=1|g" maker_opts.ctl
#fi

#Run Maker on chosen subset
export BLAST_USAGE_REPORT=FALSE
singularity exec --bind /work,/scratch,/tmp ${MAKER} \
  maker \
  -genome ${tmpDir}/genome_contigs.fasta \
  -cpus ${SLURM_CPUS_PER_TASK} \
  -TMP ${tmpDir}

sed -i "s|genome_contigs_datastore|${ARRAY_ID}/genome_contigs.maker.output/genome_contigs_datastore|g" \
  ${tmpDir}/genome_contigs.maker.output/genome_contigs_master_datastore_index.log

cp ${tmpDir}/genome_contigs.maker.output/genome_contigs_master_datastore_index.log ${outdir}/maker_array/datastore_index_${ARRAY_ID}.log
