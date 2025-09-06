# Barcode-protocol
Protocol to track single cell derived clones using DNA barcoding combined with single cell RNA sequencing

For computational analysis of barcode DNA amplicon sequencing, you will need to run the script "amplicon_analysis.sh". It requires the following files:
a. amplicon_analysis.sh
b. pair_amplicon_fastqs.py
c. Levenshtein.sh
d. Levenshtein_txt.R

To associate single cell RNA profiles from single cell RNA sequencing with DNA barcodes, you will need to follow the steps in the custom analysis pipeline "assoc_bcs_steps.sh". It requires the following files:
a. assoc_bcs_steps.sh
b. associate_barcodes_with_10x_cells.R
c. Levenshtein.R
d. fa_csv.R
