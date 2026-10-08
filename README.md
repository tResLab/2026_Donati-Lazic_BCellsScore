# B-score | Classical Hodgkin Lymphoma

We defined the B-score, a surrogate of normal B-cell content surrounding tumor cells, definining at diagnosis high-risk classical Hodgkin lymphoma (cHL) patients. 

This repository provides the computational pipeline to train the model, select the optimal hyperparameters, derived optimal coefficients and apply the derived score to individual patients in clinical practice.

## Table of Contents

- [Overview](#overview)
- [Workflow & Pipeline](#workflow--pipeline)
- [Repository structure](#repository-structure)
- [Data availability](#data-availability)
- [Requirements](#requirements)
- [Citation](#citation)

## Overview

This project implements a robust machine learning workflow (Elastic Net/Lasso) to derive a prognostic risk score starting from 18 B-cell related genes previously associated with PFS. These genes were obtained by a gene expression analysis conducted by means on nCounter PanCancer Immune Profiling Panel (Nanostring Technologies) on RNA extracted from diagnostic biopsies.

The normalized counts were scaled prior to running the model. Scaled parameters (mean and standard deviation) for each gene were saved and freezed to apply them to new patients.
The pipeline ensures high reproducibility by applying strict cross-validation fold rules during hyperparameter tuning. Once the optimal model is selected, training cohort coefficients are frozen to seamlessly extend the score calculation to new, unseen individual patients.

## Workflow & Pipeline

The pipeline consists of four main steps described in the source code:
1. **Scaled Parameters Freezing (Scaling)**.
The mean (μ) and standard deviation (sd) used to scale data for each gene are calculated and frozen. These parameters are required to correctly scale the gene expression values of any new incoming patient.
3. **Strict Cross-Validation Loop & Alpha Tuning**. 
The CV folds are strictly fixed so that every α value evaluates exactly the same patients. The pipeline extracts the minimum deviance associated with each specific α to identify the optimal model.
4. **Model generation and coefficients derivation**.
Elastic Net Regression is run with α=0.1, selected by applying the “Elbow” method. Given a comparable Cross-Validation Error to α=0.05, this configuration is associated with a lower λ value (λ=0.16), which reduces the shrinkage bias on the non-zero coefficients while selecting only 8 predictive variables instead of the 13 identified at α = 0.05.
5. **Score Computation**.
The linear score for each patient is calculated by summing the weighted gene expression values (coefficients × scaled expression). Since this is a protective score, the final result is multiplied by `-1` so that a higher score correctly aligns with a protective clinical effect (better outcome/lower risk).


## Repository structure

```text
2026_Donati-Lazic_BCellsScore/
├── README.md
├── scripts/
│   └── 0_Compute_B-cell_Score.R
└── parameters/
    ├── scaling_params_train.RData
    └── coef_vector_final.RData
```

**`0_Compute_B-cell_Score.R`** — R script for preprocessing gene expression data, selecting model parameters, fitting the Elastic Net Cox model, calculating gene coefficients, and generating the B-cell score.

The `parameters/` directory contains the trained model components required to apply the score to independent samples:

- **`scaling_params_train.RData`** — Contains `train_means` and `train_sds`, the gene-specific means and standard deviations estimated from the training cohort. These parameters are used to standardize new gene expression measurements consistently with the training data.
- **`coef_vector_final.RData`** — Contains `coef_vector`, the non-zero coefficients of the final Elastic Net Cox model, used to calculate the weighted score from the standardized expression values of the eight selected genes.


## Data availability

The input datasets used in this study are publicly available through ArrayExpress:

| Dataset | Accession number | Description |
|---|---|---|
| Training cohort | [E-MTAB-15890](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15890) | Gene expression matrix and associated clinical metadata |
| Validation cohort | [E-MTAB-15877](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15877) | Gene expression matrix and associated clinical metadata |

Clinical information, including progression-free survival (PFS) time and event status, is available in the corresponding metadata files.

The original datasets are not included in this repository.



## Requirements

The analyses were performed in R (>=4) using the following packages:

- `tidyverse`
- `survival`
- `survminer`
- `glmnet`


## Citation

If you use this code or the B-cell score in your research, please cite our paper:

- **Donati B, Lazic T, Valli R et al., " Spatial transcriptomics identifies a B-cells–depleted profile defining underlying high-risk classic Hodgkin Lymphoma.**

Full citation and DOI will be added upon publication.
