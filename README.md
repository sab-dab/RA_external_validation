# RA External Validation (Gut Microbiome Study)

## Overview
This repository contains code for external validation of a microbiome-based machine learning model for Rheumatoid Arthritis (RA).

The model was originally trained on a large 16S rRNA dataset and evaluated on an independent dataset to assess generalizability.

---

## Methods Summary
- Raw FASTQ files processed using DADA2
- Amplicon Sequence Variants (ASVs) generated
- Taxonomy assigned using SILVA (v138.1)
- Features aggregated to genus level
- Prediction performed using a pre-trained XGBoost model
- Model performance evaluated using ROC-AUC and confusion matrix

---

## Workflow
1. Download FASTQ files (external dataset)
2. Perform quality filtering and trimming
3. Learn error rates and infer ASVs
4. Merge paired reads and remove chimeras
5. Assign taxonomy (SILVA database)
6. Aggregate ASVs to genus level
7. Align features with training dataset
8. Apply trained model
9. Evaluate performance (ROC-AUC)

---

## Requirements
- R (version ≥ 4.3)
- Required packages:
  - dada2
  - dplyr
  - stringr
  - tibble
  - pROC
  - xgboost

---

## Required Files
- FASTQ files (external validation dataset)
- SILVA database (v138.1)
- Pre-trained model file (`model_xgb.model`)
- Training feature list (`train_asvs.rds`)

---

## Output
- ROC-AUC value
- Confusion matrix

---

## External Validation Dataset
Sun Y et al. (2019)  
Characteristics of gut microbiota in patients with rheumatoid arthritis in Shanghai, China  
Front Cell Infect Microbiol. 2019;9:369.

---

## Author
Sabira Dabeer

---

## Notes
This repository is provided for research and reproducibility purposes.
