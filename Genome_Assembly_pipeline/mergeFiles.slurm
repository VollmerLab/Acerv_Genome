#!/bin/bash

#SBATCH --job-name=mergeFASTQ
#SBATCH --partition=short
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G

outName=${1}; shift 1
fastqIN="$@"

cat ${fastqIN[@]} > ${outName}
