# Proteomics MI risk prediction


This repository contains code for peer review for the paper "Machine
learning-based proteomics analysis identified 23 protein biomarkers for enhanced
risk prediction of myocardial infarction in Chinese and Europeans". When this
paper is published, a link to the publication will be provided.


## Targets workflow

This analysis code uses the {targets} R package. See
<https://docs.ropensci.org/targets> for more information. In summary: `R/`
contains functions for data preparation, analysis and presentation; `Rmd/`
contains R markdown files for creating documents; `targets/` contains target
definitions; and `_targets.R` defines the overall workflow. The `run.R` script
will run all analyses.

Example analyses without using {targets} are shown in `example_analysis.R`.

## Data directories

Use environment variables RAWDATA_DIR, RAWDATA_DIR2 and UKB_RAWDATA_DIR to set
relevant data directories. Use the environment variable TARGETS_STORE to set the
targets data store directory.

## Key functions

* `prepare_participant_data()` prepares participant data for all CKB participants
from raw data
* `prepare_olink_data()`, `prepare_somascan_data()`, `scale_protein_data()`
prepare protein data and scale protein values
* `individual_proteins_analysis()` estimates associations of individual proteins
with MI
* `train_lasso_logistic()` trains a LASSO logistic regression model
* `get_predictions()` obtains predictions given a trained model and test data set
* `calculate_c_stat()`, `calculate_nri()`, `get_survfit()`, `get_calibration()`
for validation

These functions and this analysis code were developed specifically for these
analyses. They are not designed for general use.
