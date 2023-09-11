Create Phylogenetic Tree of Acroporids using Single Copy Orthologs & Orthofinder

1. Download Acropora proteomes from https://marinegenomics.oist.jp/gallery/gallery/index
```
mkdir -p /scratch/j.selwyn/phylogenetics/slurm
sbatch \
  --output=/scratch/j.selwyn/phylogenetics/slurm/downloadFiles_%j.output \
  /work/vollmer/software/jds_scripts/runRscript.slurm \
    /work/vollmer/software/jds_scripts/r_utils/download_acroporid_proteomes.R \
    /scratch/j.selwyn/phylogenetics/raw_data
```
2. Download Other proteomes used for annotation
```
srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
source activate ncbi_datasets
storeDir=/work/vollmer/Databases/reference_proteomes
mkdir -p ${storeDir}
cd ${storeDir}
finalDir=/scratch/j.selwyn/phylogenetics/raw_data

#Stylophora pistillata
datasets download genome accession GCF_002571385.1 --include gff3,genome --filename GCF_002571385.1.zip
unzip GCF_002571385.1.zip
mv ncbi_dataset/data/GCF_002571385.1 ./
rm -rf ncbi_dataset; rm -rf GCF_002571385.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_002571385.1/*.fna >\
  ${finalDir}/spis.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_002571385.1/genomic.gff >\
  ${finalDir}/spis.gff.gz


#Pocillopora damicornis
datasets download genome accession GCF_003704095.1 --include gff3,genome --filename GCF_003704095.1.zip
unzip GCF_003704095.1.zip
mv ncbi_dataset/data/GCF_003704095.1 ./
rm -rf ncbi_dataset; rm -rf GCF_003704095.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_003704095.1/*.fna >\
  ${finalDir}/pdam.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_003704095.1/genomic.gff >\
  ${finalDir}/pdam.gff.gz


#Actinia tenebrosa
datasets download genome accession GCF_009602425.1 --include gff3,genome --filename GCF_009602425.1.zip
unzip GCF_009602425.1.zip
mv ncbi_dataset/data/GCF_009602425.1 ./
rm -rf ncbi_dataset; rm -rf GCF_009602425.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_009602425.1/*.fna >\
  ${finalDir}/aten.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_009602425.1/genomic.gff >\
  ${finalDir}/aten.gff.gz


#Nematostella vectensis
datasets download genome accession GCF_000209225.1 --include gff3,genome --filename GCF_000209225.1.zip
unzip GCF_000209225.1.zip
mv ncbi_dataset/data/GCF_000209225.1 ./
rm -rf ncbi_dataset; rm -rf GCF_000209225.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_000209225.1/*.fna >\
  ${finalDir}/nvec.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_000209225.1/genomic.gff >\
  ${finalDir}/nvec.gff.gz


#Exaiptasia diaphana
datasets download genome accession GCF_001417965.1 --include gff3,genome --filename GCF_001417965.1.zip
unzip GCF_001417965.1.zip
mv ncbi_dataset/data/GCF_001417965.1 ./
rm -rf ncbi_dataset; rm -rf GCF_001417965.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_001417965.1/*.fna >\
  ${finalDir}/edia.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_001417965.1/genomic.gff >\
  ${finalDir}/edia.gff.gz
```
3. Download last batch of proteomes from NCBI used in hexacoral Genomes
```
source activate ncbi_datasets
storeDir=/work/vollmer/Databases/reference_proteomes
mkdir -p ${storeDir}
cd ${storeDir}
finalDir=/scratch/j.selwyn/phylogenetics/raw_data

#A digitifera
datasets download genome accession GCF_000222465.1 --include gff3,genome --filename GCF_000222465.1.zip
unzip GCF_000222465.1.zip
mv ncbi_dataset/data/GCF_000222465.1 ./
rm -rf ncbi_dataset; rm -rf GCF_000222465.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_000222465.1/*.fna >\
  ${finalDir}/adig_shinzato1.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_000222465.1/genomic.gff >\
  ${finalDir}/adig_shinzato1.gff.gz

#A hyacithus
datasets download genome accession GCA_020536085.1 --include gff3,genome --filename GCA_020536085.1.zip
unzip GCA_020536085.1.zip
mv ncbi_dataset/data/GCA_020536085.1 ./
rm -rf ncbi_dataset; rm -rf GCA_020536085.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_020536085.1/*.fna >\
  ${finalDir}/ahya_palumbi.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_020536085.1/genomic.gff >\
  ${finalDir}/ahya_palumbi.gff.gz

#A millepora
datasets download genome accession GCF_004143615.1 --include gff3,genome --filename GCF_004143615.1.zip
unzip GCF_004143615.1.zip
mv ncbi_dataset/data/GCF_004143615.1 ./
rm -rf ncbi_dataset; rm -rf GCF_004143615.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_004143615.1/*.fna >\
  ${finalDir}/amil_ying.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_004143615.1/genomic.gff >\
  ${finalDir}/amil_ying.gff.gz

#A millepora
datasets download genome accession GCF_013753865.1 --include gff3,genome --filename GCF_013753865.1.zip
unzip GCF_013753865.1.zip
mv ncbi_dataset/data/GCF_013753865.1 ./
rm -rf ncbi_dataset; rm -rf GCF_013753865.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_013753865.1/*.fna >\
  ${finalDir}/amil_fuller.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCF_013753865.1/genomic.gff >\
  ${finalDir}/amil_fuller.gff.gz

#M capitata
datasets download genome accession GCA_006542545.1 --include gff3,genome --filename GCA_006542545.1.zip
unzip GCA_006542545.1.zip
mv ncbi_dataset/data/GCA_006542545.1 ./
rm -rf ncbi_dataset; rm -rf GCA_006542545.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_006542545.1/*.fna >\
  ${finalDir}/mcap.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_006542545.1/genomic.gff >\
  ${finalDir}/mcap.gff.gz

#M sp.
datasets download genome accession GCA_023091755.1 --include gff3,genome --filename GCA_023091755.1.zip
unzip GCA_023091755.1.zip
mv ncbi_dataset/data/GCA_023091755.1 ./
rm -rf ncbi_dataset; rm -rf GCA_023091755.1.zip; rm -rf README.md
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_023091755.1/*.fna >\
  ${finalDir}/msp.fasta.gz
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  ${storeDir}/GCA_023091755.1/genomic.gff >\
  ${finalDir}/msp_fuller.gff.gz
```
4. Put *A. cervicornis* into raw data
```
pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  /work/vollmer/k2_nanopore_genome/k2_structuralAnnotations.gff3 >\
  /scratch/j.selwyn/phylogenetics/raw_data/acer.gff.gz

pigz -p${SLURM_CPUS_PER_TASK} -k -c \
  /work/vollmer/k2_nanopore_genome/k2_genome.fasta >\
  /scratch/j.selwyn/phylogenetics/raw_data/acer.fasta.gz
```
5. Go through and delete duplicate species & species without gff files
```
DIR=/scratch/j.selwyn/phylogenetics
rm \
  ${DIR}/raw_data/adig_shinzato1* \
  ${DIR}/raw_data/amil_ying* \
  ${DIR}/raw_data/ahya_palumbi* \
  ${DIR}/raw_data/mcap* \
  ${DIR}/raw_data/msp*
```
6. Extract Proteins from FASTA/GFF files
```
sbatch \
  --dependency=afterany:38171237 \
  --output=/scratch/j.selwyn/phylogenetics/slurm/proteinExtract_%j.output \
  /work/vollmer/software/jds_scripts/extractProteinGFF.slurm \
    /scratch/j.selwyn/phylogenetics/proteins \
    /scratch/j.selwyn/phylogenetics/raw_data

#When finished run
grep -c '^>' /scratch/j.selwyn/phylogenetics/proteins/*.fasta > /scratch/j.selwyn/phylogenetics/proteins/protein_counts.txt
```
7. Run BUSCO on all Proteins
```
sbatch \
  --dependency=afterany:38171349 \
  --output=/scratch/j.selwyn/phylogenetics/slurm/buscoRun_%j.output \
  /work/vollmer/software/jds_scripts/get_busco_proteins.slurm \
    /scratch/j.selwyn/phylogenetics/buscos \
    /scratch/j.selwyn/phylogenetics/proteins
```
8. Run Orthofinder (will automatically run second orthofinder for MSA/FastTree & TimeCalibration/CAFE)
```
sinfo -p short --Format=time,nodes,cpus,socketcorethread,memory,nodeai,features

sbatch \
  --cpus-per-task=28 \
  --mem=250G \
  --output=/scratch/j.selwyn/phylogenetics/slurm/OrthoFinder_%j.output \
  /work/vollmer/software/jds_scripts/orthofinder.slurm \
    /scratch/j.selwyn/phylogenetics/orthofinder \
    /scratch/j.selwyn/phylogenetics/proteins
```

9. Get Protein Annotations
```
srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
outdir=/scratch/j.selwyn/phylogenetics
mkdir -p ${outdir}/annotations/slurm
cd ${outdir}/annotations

# 1. merge all protein files into one mega fasta
cat ${outdir}/proteins/*.fasta > ${outdir}/annotations/combined_proteins.fasta

# 2. remove duplicated gene names
awk '/^>/{f=!d[$1];d[$1]=1}f' ${outdir}/annotations/combined_proteins.fasta >\
  ${outdir}/annotations/combined_proteins_dedupe.fasta

# 3. run SWISSPROT/BLAST on all genes
sbatch \
   --job-name=BLAST \
   --output=${outdir}/annotations/slurm/blast_%A_%a.output \
   --array=0-$((250-1))%40 \
   /work/vollmer/software/jds_scripts/functionalAnnotateBLAST.slurm \
   ${outdir}/annotations \
   ${outdir}/annotations/combined_proteins_dedupe.fasta \
   250 \
   10

# 4. gather array results into useable format - somthing goes wrong with parsing when simply cat together
sbatch \
  --dependency=afterany:38176209 \
  --output=/scratch/j.selwyn/phylogenetics/annotations/slurm/blastReformat_%j.output \
  /work/vollmer/software/jds_scripts/runRscript.slurm \
    /work/vollmer/software/jds_scripts/r_utils/reformat_blast_kegg.R \
    /scratch/j.selwyn/phylogenetics/annotations/blast_array
```
10. Join Orthofinder & Annotation Results
```
sbatch \
  --output=/scratch/j.selwyn/phylogenetics/slurm/ortho_annotate_%j.output \
  /work/vollmer/software/jds_scripts/runRscript.slurm \
  /work/vollmer/software/jds_scripts/r_utils/annotate_orthofinder.R \
    /scratch/j.selwyn/phylogenetics/annotations/kegg_annotations.csv.gz \
    /scratch/j.selwyn/phylogenetics/orthofinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv \
    /scratch/j.selwyn/phylogenetics
```

11. Time Calibration & Cafe5
- Put fossil calibration file in folder `/scratch/j.selwyn/phylogenetics/fossil_estimates.txt`
```
sbatch \
  --output=/scratch/j.selwyn/phylogenetics/slurm/time_calibration_cafe_%j.output \
  /work/vollmer/software/jds_scripts/runCAFE.slurm \
    /scratch/j.selwyn/phylogenetics/time_calibration_cafe \
    /scratch/j.selwyn/phylogenetics/orthofinder/*/Species_Tree/SpeciesTree_rooted.txt \
    /scratch/j.selwyn/phylogenetics/fossil_estimates.txt \
    /scratch/j.selwyn/phylogenetics/orthofinder/*/Orthogroup_Sequences \
    /scratch/j.selwyn/phylogenetics/orthofinder_annotations.csv.gz \
    100
```
