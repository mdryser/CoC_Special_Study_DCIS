# Code Repository for BMJ Study  
**Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS)**  
*Marc D. Ryser, E. Shelley Hwang, et al. BMJ, 2025*

## Overview

This repository contains R code for the analyses published in:

> **Ryser MD, Hwang ES, et al.**  
> *Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS).*  
> BMJ, 2025.

The code generates figures and tables describing absolute and relative risks of ipsilateral invasive breast cancer (iIBC) and disease-specific survival (DSS) among women managed without immediate surgery for DCIS.

> ⚠️ **Data Note**: The analytic dataset (`analytic_cohort.RData`) and life table input are not publicly available. This repository provides code only, for transparency and reproducibility.

---

## Repository Contents

### `BMJ_2025_Ryser_Hwang_Absolute_Risks.R`

Generates primary results on absolute risk:
- **Table 1**: Patient and tumor characteristics
- **Figures 2A–B**: Cumulative incidence of iIBC
- **Figures 3A–B**: Disease-specific survival
- **Supplementary Figures A & C**: Time to surgery and competing risks

### `BMJ_2025_Ryser_Hwang_Relative_Risks.R`

Generates relative risk estimates and sensitivity analyses:
- **Supplementary Figure B**: Partial residual plot (age vs. iIBC risk)
- **Supplementary Table E**:
  - A: Univariate Cox models (baseline covariates)
  - B: Time-dependent Cox models with fixed ET durations
  - C: Monte Carlo–sampled ET durations in time-dependent Cox models

### `aux.fun.R`

Provides supporting functions:
- `Risk.KM()` – Absolute risk estimation from Kaplan-Meier fits  
- `DSS.KM()` – Disease-specific survival estimation  
- `OC.surv()` – Expected mortality modeling using U.S. life tables

---

## Software Requirements

The following R packages are required:

```r
install.packages(c(
  "haven", "dplyr", "survival", "ggplot2", "ggsci", "cmprsk",
  "survminer", "ggsurvfit", "tidycmprsk", "gtsummary", "flextable",
  "mice", "mitools"
))
