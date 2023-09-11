#!/bin/bash

#SBATCH --job-name=postAssembly
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outdir=${1}
draft_genome=${2}
genome_name=${3}
min_length=${4}

scriptDir=/work/vollmer/software/jds_scripts

mkdir -p ${outdir}
cd ${outdir}

#Move Mitochondrial Genome to main directory
mv ${outdir}/mitogenome/*/mt_assembly/*megahit.result/*mitogenome.fa ${outdir}/${genome_name}_mitogenome.fasta
mv ${outdir}/mitogenome/*/mt_annotation/*fa.result/*.gbf ${outdir}/${genome_name}_mitogenome.gbf
mv ${outdir}/mitogenome/*/mt_annotation/*fa.result/circos.png ${outdir}/${genome_name}_mitogenome.png
mv ${outdir}/mitogenome/*/mt_annotation/*fa.result/circos.svg ${outdir}/${genome_name}_mitogenome.svg

# module load anaconda3; source activate bioperl
# bp_genbank2gff3.pl

#Estimate Genome Size
module load R/4.0.3
genome_size=$(Rscript \
  ${scriptDir}/r_utils/jellyfish_genome_size.R \
    ${outdir}/jellyfish/*.histo)
echo ${genome_size}

#Sort genome contigs by length then rename scaffolds then filter those shorter than ${min_length}
module load anaconda3; source activate bioawk
bioawk -v LENGTH="${min_length}" -c fastx '{ if(length($seq) > LENGTH) { print ">"$name; print $seq }}' \
  ${draft_genome} > ${outdir}/tmp_lengthFilter.fasta

conda deactivate; module load anaconda3; source activate bbmap
sortbyname.sh in=${outdir}/tmp_lengthFilter.fasta out=${outdir}/tmp_sort.fasta length descending
rename.sh in=${outdir}/tmp_sort.fasta out=${outdir}/tmp_sort_rename.fasta prefix=Acerv_scaffold

mv ${outdir}/tmp_sort_rename.fasta ${outdir}/${genome_name}.${draft_genome##*.}
rm ${outdir}/tmp_*fasta

# Make Table with contig lengths
conda deactivate; module load anaconda3; source activate bioawk
bioawk -c fastx '{ print $name, length($seq) }' < ${outdir}/${genome_name}.${draft_genome##*.} > ${outdir}/${genome_name}_scaffoldLengths.tsv

# Run assembly assessment on final assembly
conda deactivate
bash \
  ${scriptDir}/assessAssembly.slurm \
    ${outdir}/${genome_name}.${draft_genome##*.}

# Make json for interactive plot
#After making json transfer to forked: https://github.com/ammaraziz/assembly-stats and follow steps there to produce github page with summary stats/plots
module load anaconda3; source activate plot_assembly_stats
asm2stats.minmaxgc.pl ${outdir}/${genome_name}.${draft_genome##*.} > ${outdir}/output.assembly-stats.json
conda deactivate

# Make Final Plots
Rscript \
  ${scriptDir}/r_utils/summarize_post_assembly.R \
    ${outdir} \
    ${genome_size}
