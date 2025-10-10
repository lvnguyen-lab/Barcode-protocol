#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --mem=180G
#SBATCH -J assoc_bc_steps
#SBATCH -p veryhimem
#SBATCH -c 10
#SBATCH -N 1
#SBATCH --mail-type=END
#SBATCH --mail-user=sarah.kronheim@uhn.ca
#SBATCH -o %x-%j.out

# Load specific run's parameters 
WD=$(pwd -P)
BASE_DIR=/cluster/projects/nguyengroup/scRNAseq
RUN_NAME=$1

if [ "$WD" != "$BASE_DIR/runs/$RUN_NAME" ]
then
        echo "Error: you should run this from the specific run directory ($BASE_DIR/runs/$RUN_NAME)"
        exit 1
fi

META_FILE=$(ls -1 $BASE_DIR/runs/$RUN_NAME/metadata/*csv)
SAMPLE_NAME=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $4 }' | tr -d '"')
SAMPLE_ID=$(tail -n +2 $META_FILE | sed -n "$SLURM_ARRAY_TASK_ID"p | awk -F',' '{ print $2 }' | tr -d '"')

mkdir -p ${WD}/barcode_assoc/${SAMPLE_NAME}/BWA_index

ASSOC_DIR=${WD}/barcode_assoc/${SAMPLE_NAME}

##Associate GFP_barcode with cellid from 10X sequencing
# This script calls both Levenshtein.R and associate_barcodes_with_10x_cells.R scripts
# Make variables for PATHs and for sample names
    # Sample names -> can use SAMPLE_ID(=SITTXXX) from cellranger input
# If add to pipeline, will only be 1 sample at a time so no need of loops

# Pool reads for individual sample
cat ${WD}/fastq/raw/${SAMPLE_ID}_*_R2_001.fastq.gz > ${WD}/fastq/combined/${SAMPLE_ID}_direct.s_1.r_2.fq.gz; 
rm ${WD}/fastq/raw/${SAMPLE_ID}_*

#Use Cudadapt to trim Read2 file based on constant up/downstream sequences, and filter based on size (31bp) & quality (Phred ≥30)
cutadapt -a TCGTTTTACA...CTGTCGACCT --discard-untrimmed -q 30 -m 31 -M 31 -o ${ASSOC_DIR}/${SAMPLE_ID}_direct_r2.fq ${WD}/fastq/combined/${SAMPLE_ID}_direct.s_1.r_2.fq.gz
rm ${WD}/fastq/combined/${SAMPLE_ID}_*

#Use grep to filter fastq file based on expected barcode sequence pattern 
grep -A 2 -B 1 '..GTT..ACC..TTT..ATC..GAT..CAGT\|..GTT..ACC..TTT..ATC..GAT..GTAC' ${ASSOC_DIR}/${SAMPLE_ID}_direct_r2.fq | sed '/--/d' > ${ASSOC_DIR}/${SAMPLE_ID}_direct_filtered.fq
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_r2.fq

#Get reverse complement fastq file
fastx_reverse_complement -i ${ASSOC_DIR}/${SAMPLE_ID}_direct_filtered.fq -o ${ASSOC_DIR}/${SAMPLE_ID}_direct_rc.fastq
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_filtered.fq

#Convert fastq --> seqs --> csv files
cat ${ASSOC_DIR}/${SAMPLE_ID}_direct_rc.fastq | sed -n '2 ~ 4p' > ${ASSOC_DIR}/${SAMPLE_ID}_direct.seqs

cat ${ASSOC_DIR}/${SAMPLE_ID}_direct.seqs > ${ASSOC_DIR}/${SAMPLE_ID}_direct.csv
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct.seqs

#Use R to tally unique read counts from csv files
Rscript ${BASE_DIR}/scripts/Levenshtein.R ${ASSOC_DIR}/${SAMPLE_ID}_direct.csv ${ASSOC_DIR}/${SAMPLE_ID}_direct_BC1.csv
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct.csv

# Add step - convert this output csv to fasta format - this will go into bwa 
Rscript ${BASE_DIR}/scripts/fa_csv.R ${ASSOC_DIR} ${SAMPLE_ID}_direct_BC1.csv ${SAMPLE_ID}_direct_BC1.fa
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_BC1.csv

# Makes reference file - align original fastq to reference - need info abt cell barcode index

# Generate BWA alignment index files
cd ./barcode_assoc/${SAMPLE_NAME}/BWA_index

bwa index -p BC1_${SAMPLE_ID}_didx -a bwtsw ../${SAMPLE_ID}_direct_BC1.fa

#BWA alignment - run in WD
cd ${WD}

bwa aln -n 0 -o 0 -e 0 -l 10 -k 0 -t 4 ${ASSOC_DIR}/BWA_index/BC1_${SAMPLE_ID}_didx ${ASSOC_DIR}/${SAMPLE_ID}_direct_rc.fastq > ${ASSOC_DIR}/${SAMPLE_ID}_direct.sai

#Convert .sai to .sam
bwa samse ${ASSOC_DIR}/BWA_index/BC1_${SAMPLE_ID}_didx ${ASSOC_DIR}/${SAMPLE_ID}_direct.sai ${ASSOC_DIR}/${SAMPLE_ID}_direct_rc.fastq > ${ASSOC_DIR}/${SAMPLE_ID}_direct.sam
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_rc.fastq
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct.sai

#Convert .sam to .bam
samtools view -b ${ASSOC_DIR}/${SAMPLE_ID}_direct.sam -o ${ASSOC_DIR}/${SAMPLE_ID}_direct.bam

#Filter out only mapped reads
#-f 4 will give you unmapped reads output
#-F 4 will give you mapped reads output
#-b is for bam file output, otherwise default is as sam file

#Keep only mapped reads, and convert BAM to SAM
samtools view -b -F 4 ${ASSOC_DIR}/${SAMPLE_ID}_direct.bam > ${ASSOC_DIR}/${SAMPLE_ID}_direct_mapped.bam
samtools view -F 4 ${ASSOC_DIR}/${SAMPLE_ID}_direct.bam > ${ASSOC_DIR}/${SAMPLE_ID}_direct_mapped.sam
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct.bam

#To see number of mapped reads
samtools flagstat ${ASSOC_DIR}/${SAMPLE_ID}_direct.sam
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct.sam

#Convert file names to number
mv ${ASSOC_DIR}/${SAMPLE_ID}_direct_mapped.sam ${ASSOC_DIR}/${SAMPLE_NAME}_direct_mapped.sam

#Run published custom script against 10X data to get a table of cell BC with internal GFP BC
#Run from /mnt/scratchb/cclab/nguyen03/ directory
# 3 names each - tumor, SITT (index), S1-36

Rscript ${BASE_DIR}/scripts/associate_barcodes_with_10x_cells.R ${SAMPLE_NAME} ${WD} ${ASSOC_DIR}/${SAMPLE_NAME}_direct_mapped.sam ${WD}/report


# Remove remaining intermediate files, keeping only csv files
rm -r ${ASSOC_DIR}/BWA_index
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_BC1.fa
rm ${ASSOC_DIR}/${SAMPLE_ID}_direct_mapped.bam
rm ${ASSOC_DIR}/${SAMPLE_NAME}_direct_mapped.sam
rm ${WD}/fastq/raw/${SAMPLE_ID}*