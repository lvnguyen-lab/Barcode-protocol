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
#   - Operator folder and script structure:
#       - the operator has set up a folder for amplicon analysis, in which this shell script is present
#       - this folder has a folder for "runs" and a folder for "scripts", and these folders are named as such
#       - the "scripts" folder contains the "pair_amplicon_fastqs.py" python script and the "Levenshtein_txt.R" R script
# Outputs
#   - csv files with sequences and reads for pre- and post-Levenshtein aggregation for each sample and each library
#     - Be aware that all 3 library sequences are tested against each sample, so it is critical to know exactly which library any given sample should belong to

# Setup
BASE_DIR=/cluster/projects/nguyengroup/
WD=$(pwd -P)
AMP_RUN=$1
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
  echo "Please enter the type of reads you have want to perform"
  echo "Enter 'paired' for paired reads"
  echo "Enter 'unpaired' for unpaired reads"
  exit 1
fi

if [[ $3 == "" ]];
  then
    QUALITY=30
fi

RUN_DIR=${WD}/runs/${AMP_RUN} 
AMP_DIR=${RUN_DIR}/fastq 

mkdir -p ${RUN_DIR}/Cutadapt
mkdir -p ${RUN_DIR}/Grep
mkdir -p ${RUN_DIR}/Final_counts
mkdir -p ${RUN_DIR}/Summary_stats

# Activate environment
conda init zsh
source ~/.zshrc
conda activate NMF_env

# Run Python script to pair R1 and R2 files in the AMP_DIR; saves a CSV file with columns for R1 and R2 that then must be read
python3 ${WD}/scripts/pair_amplicon_fastqs.py ${AMP_DIR} ${RUN_DIR}

# Load in paired csv file
PAIRED_FASTQS=${RUN_DIR}/paired_fastq.csv

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
    cutadapt -a TAGAAGGCACAGGTCGACAG...TGTTCTAGACATATGGTCTCTCCTAG -A AGGAGAGACCATATGTCTAGAACA...CTGTCGACCTGTGCCTT --discard-untrimmed -o ${RUN_DIR}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq -p ${RUN_DIR}/Cutadapt/${sample}_R2_trimmed_q${QUALITY}_${TYPE}.fq ${AMP_DIR}/${R1} ${AMP_DIR}/${R2} > ${RUN_DIR}/Summary_stats/${sample}_adapter_trim_${TYPE}.txt

    cutadapt -m 40 -M 40 -q ${QUALITY} -o ${RUN_DIR}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq -p ${RUN_DIR}/Cutadapt/${sample}_rev_q${QUALITY}_${TYPE}.fq ${RUN_DIR}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq ${RUN_DIR}/Cutadapt/${sample}_R2_trimmed_q${QUALITY}_${TYPE}.fq > ${RUN_DIR}/Summary_stats/${sample}_quality_trim_${TYPE}.txt
  else
    cutadapt -a TAGAAGGCACAGGTCGACAG...TGTTCTAGACATATGGTCTCTCCTAG --discard-untrimmed -o ${RUN_DIR}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq ${AMP_DIR}/${R1} > ${RUN_DIR}/Summary_stats/${sample}_adapter_trim_${TYPE}.txt

    cutadapt -m 40 -M 40 -q ${QUALITY} -o ${RUN_DIR}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq ${RUN_DIR}/Cutadapt/${sample}_R1_trimmed_q${QUALITY}_${TYPE}.fq > ${RUN_DIR}/Summary_stats/${sample}_quality_trim_${TYPE}.txt
  fi

  #Run grep to retrieve the reads that match the expected sequence

  grep -A 2 -B 1 'ATCG..GAT..ATC..GGT..ATA..TGT..AAC..ATCG' ${RUN_DIR}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq | sed -n '/^ATCG/p' > ${RUN_DIR}/Grep/${sample}_BC1_q${QUALITY}_${TYPE}.txt

  grep -A 2 -B 1 'GACC..GAT..ATC..GGT..ATA..TGT..AAC..GACC' ${RUN_DIR}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq | sed -n '/^GACC/p' > ${RUN_DIR}/Grep/${sample}_BC2_q${QUALITY}_${TYPE}.txt

  grep -A 2 -B 1 'TGGT..CTA..TAG..CGA..AAA..TAG..ATC..TGGT' ${RUN_DIR}/Cutadapt/${sample}_fwd_q${QUALITY}_${TYPE}.fq | sed -n '/^TGGT/p' > ${RUN_DIR}/Grep/${sample}_BC3_q${QUALITY}_${TYPE}.txt

  # Run Levenshtein.R to group single mismatches
  echo "Running Levenshtein script on BC1"
  sbatch Levenshtein.sh ${RUN_DIR}/Grep/${sample}_BC1_q${QUALITY}_${TYPE}.txt ${RUN_DIR}/Final_counts/${sample}_BC1_q${QUALITY}_${TYPE}

  echo "Running Levenshtein script on BC2"
  sbatch Levenshtein.sh ${RUN_DIR}/Grep/${sample}_BC2_q${QUALITY}_${TYPE}.txt ${RUN_DIR}/Final_counts/${sample}_BC2_q${QUALITY}_${TYPE}

  echo "Running Levenshtein script on BC3"
  sbatch Levenshtein.sh ${RUN_DIR}/Grep/${sample}_BC3_q${QUALITY}_${TYPE}.txt ${RUN_DIR}/Final_counts/${sample}_BC3_q${QUALITY}_${TYPE}

done < ${PAIRED_FASTQS}

