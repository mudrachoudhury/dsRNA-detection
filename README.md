# dsRNA_detection
Detect dsRNAs from known RNA Editing sites in REDI Portal (2020). Written to run on SLURM HPC

Used in the following publication:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03012-w

written by: Mudra Choudhury

Table of contents:

STEP 1: Format all REDI portal sites

STEP 2: Identify Editing Enriched Regions (EERs)

STEP 3: Merge overlapping/nearby EERs (1kb used)

STEP 4: Apply minimum and maximum length cutoffs to merged EERs

STEP 5: Run RNAfold to retrieve Minimum Free Energy (MFE) and base pairing

STEP 6: Apply stem length (double stranded length) and mismatch cutoff




