#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --mem=100G
#SBATCH -J sk_amplicon_processing
#SBATCH -p veryhimem
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --mail-type=END
#SBATCH --mail-user=sarah.kronheim@uhn.ca
#SBATCH -o %x-%j.out

# USAGE ASSUMPTIONS: 
#   - fastq files:
#       - fastq data is in the ampliconData folder in the nguyengroup 
#       - fastq data is set up with run_name/fastq and run_name/QC folders
#       - the run_name is used with this shell script
#           - when this script is called, type "sbatch amplicon_seq_processing_adapted.sh run_name" into the terminal
#       - ONLY fastq files are in the fastq folder
#       - data files are in .fastq.gz format
#       - The reads are paired, so R1 and R2 files are present for each sample 
#   - Operator folder and script structure:
#       - the operator has set up a folder for amplicon analysis, in which this shell script is present
#       - this folder has a folder for "runs" and a folder for "scripts", and these folders are named as such
#       - the "scripts" folder contains the "pair_amplicon_fastqs.py" python script and the "Levenshtein_txt.R" R script
# Outputs
#   - csv files with sequences and reads for pre- and post-Levenshtein aggregation for each sample and each library
#     - Be aware that all 3 library sequences are tested against each sample, so it is critical to know exactly which library any given sample should belong to

# Setup
BASE_DIR=/cluster/projects/nguyengroup/scRNAseq
WD=$(pwd -P)
RUN_NAME=$1
TYPE=$2
QUALITY=$3

if [[ $1 == "" ]];
  then
    echo "Error: missing run name"
    exit 1
fi

if [[ $2 == "paired" ]] || [[ $2 == "unpaired" ]];
then 
  TYPE=$2
  echo ${TYPE}
else
  echo "Please enter whether you want to analyze paired reads together or as unpaired"
  echo "Enter 'paired' for paired reads"
  echo "Enter 'unpaired' for unpaired reads"
  exit 1
fi

if [[ $3 == "" ]];
  then
    QUALITY=30
fi

if [ "$WD" != "$BASE_DIR/runs/$RUN_NAME" ]
then
        echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/$RUN_NAME)"
        exit 1
fi

AMP_DIR=${WD}/fastq/amplicon 

mkdir -p ${WD}/Cutadapt
mkdir -p ${WD}/Grep
mkdir -p ${WD}/Final_counts
mkdir -p ${WD}/Summary_stats

# Activate environment
conda init zsh
source ~/.zshrc
conda activate NMF_env

# Run Python script to pair R1 and R2 files in the AMP_DIR; saves a CSV file with columns for R1 and R2 that then must be read
python3 ${BASE_DIR}/scripts/pair_amplicon_fastqs.py ${AMP_DIR} ${WD}

# Load in paired csv file
PAIRED_FASTQS=${WD}/paired_fastq.csv

# Load modules
module load R
module load cutadapt

# Read csv and run script on paired reads
while IFS=, read -r sample R1 R2
do  
  echo "analyzing sample ${sample}"  

  if [[ ${TYPE} == "paired" ]]
  then
  #Run Cutadapt to trim constant regions up and downstream of barcode sequence, and also filter based on size (40bp) and quality (Phred ≥30)
  # Cut and keep barcode sequences only (this sequence works for all staggered primers, and has been adapted to the new barcode libraries made by May and naming conventions of PMCRT seq center) 
    cutadapt -a TAGAAGGCACAGGTCGACAG...TGTAAAACGACGGCCAGTG -A GACTCACTGGCCGTCGTTTTACA...CTGTCGACCTGTGCCT --discard-untrimmed -o ${WD}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq -p ${WD}/Cutadapt/${sample}_R2_trimmed_q${QUALITY}_${TYPE}.fq ${AMP_DIR}/${R1} ${AMP_DIR}/${R2} > ${WD}/Summary_stats/${sample}_adapter_trim_${TYPE}.txt

    cutadapt -m 40 -M 40 -q ${QUALITY} -o ${WD}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq -p ${WD}/Cutadapt/${sample}_rev_q${QUALITY}_${TYPE}.fq ${WD}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq ${WD}/Cutadapt/${sample}_R2_trimmed_q${QUALITY}_${TYPE}.fq > ${WD}/Summary_stats/${sample}_quality_trim_${TYPE}.txt
  else
    cutadapt -a TAGAAGGCACAGGTCGACAG...TGTAAAACGACGGCCAGTG --discard-untrimmed -o ${WD}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq ${AMP_DIR}/${R1} > ${WD}/Summary_stats/${sample}_adapter_trim_${TYPE}.txt

    cutadapt -m 40 -M 40 -q ${QUALITY} -o ${WD}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq ${WD}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq > ${WD}/Summary_stats/${sample}_quality_trim_${TYPE}.txt
  fi

  #Run grep to retrieve the reads that match the expected sequence
  grep -A 2 -B 1 'ACTG..ATC..GAT..AAA..GGT..AAC..' ${WD}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq | sed -n '/^ATCG/p' > ${WD}/Grep/${sample}_BC1_q${QUALITY}_${TYPE}.txt

  # Run Levenshtein.R to group single mismatches
  echo "Running Levenshtein script"
  sbatch ${BASE_DIR}/scripts/Levenshtein.sh ${WD}/Grep/${sample}_BC1_q${QUALITY}_${TYPE}.txt ${WD}/Final_counts/${sample}_BC1_q${QUALITY}_${TYPE}

done < ${PAIRED_FASTQS}

