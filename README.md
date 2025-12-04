# Barcode-protocol
As published in Shin HJ et al. Protocol to track single-cell-derived clones using DNA barcoding combined with single-cell RNA sequencing. STAR Protocols 2025. PMID: 41317321, DOI: 10.1016/j.xpro.2025.104229

For computational analysis of barcode DNA amplicon sequencing, you will need to run the script "amplicon_analysis.sh". It requires the following files:
a. amplicon_analysis.sh
b. pair_amplicon_fastqs.py
c. Levenshtein.sh
d. Levenshtein_txt.R

To associate single cell RNA profiles from single cell RNA sequencing with DNA barcodes, you will need to follow the steps in the custom analysis pipeline "assoc_bcs_steps.sh". It requires the following files:
a. associate_barcodes_with_10x_cells.R
b. Levenshtein.R
c. fa_csv.R

The file SLX-22131.HL2NCDSX3.s_1.contents.csv is an example metadata file. There should be only one metadata file in location "$BASE_DIR/runs/$RUN_NAME/metadata". All samples for the run can be included in this one file.
