# AAFC_Bjuncea_Salinity
R scripts and datasets supporting: “From germination to maturity: gene-bank enabled physiological screening of salinity tolerance in Brassica juncea.”

All analyses are fully reproducible by running the scripts in numerical order.

Repository Structure
AAFC_Bjuncea_Salinity/
│
├── data/
│   ├── master_run_1.xlsx
│   ├── master_run_2.xlsx
│   └── master_agronomy.xlsx
│
├── scripts/
│   ├── 01_photosynq_analysis.R
│   ├── 02_leafspec_analysis.R
│   └── 03_composite_scores_and_sti.R
│
├── outputs/
│   └── (generated automatically when scripts are run)
│
└── AAFC_Bjuncea_Salinity.Rproj

How to Reproduce Analyses

Open the R project file (AAFC_Bjuncea_Salinity.Rproj)

Run scripts in order:

source("scripts/01_photosynq_analysis.R")
source("scripts/02_leafspec_analysis.R")
source("scripts/03_combined_analysis.R")

Outputs will be written automatically to the outputs/ folder.

Key Outputs

PhotosynQ t-tests and heritability estimates

LeafSpec t-tests and heritability estimates

Composite scores (weighted by developmental stage)

Stress Tolerance Index (STI)

Composite Scores

Software Requirements

R (≥ 4.2 recommended)

Packages:

readxl

dplyr

tidyr

purrr

tibble

lme4

emmeans

FactoMineR

missMDA

ggplot2
