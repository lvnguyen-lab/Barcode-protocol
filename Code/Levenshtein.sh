#!/bin/bash
#SBATCH -t 21-00:00:00
#SBATCH --mem=5G
#SBATCH -J sk_levenshtein_amplicon_processing
#SBATCH -p long
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --mail-type=END
#SBATCH --mail-user=sarah.kronheim@uhn.ca
#SBATCH -o %x-%j.out

# Setup
INPUT=$1
OUTPUT=$2
WD=$(pwd -P)

# Load modules
module load R

# Run Levenshtein on input file
Rscript ${WD}/scripts/Levenshtein_txt.R ${INPUT} ${OUTPUT}