# Code Repository for BMJ Study  
**Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS)**  
*Marc D. Ryser, ..., E. Shelley Hwang BMJ, 2025*

## Overview

This repository contains R code for the analyses published in:

> **Ryser MD,..., Hwang ES**  
> *Cancer outcomes in women without upfront surgery for ductal carcinoma in situ (DCIS).*  
> BMJ, 2025.

The code generates figures and tables describing absolute and relative risks of ipsilateral invasive breast cancer (iIBC) and disease-specific survival (DSS) among women managed without immediate surgery for DCIS.

> ⚠️ **Data Note**: The analytic dataset and life table input required to run these scripts are **not included** in this repository.

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
```

---

## 🔍 Input Files (Not Included)

To run the code, the following input files must be available in the working directory. These are **not distributed publicly** due to data use restrictions.

### 1. `analytic_cohort.RData`

This file contains the pre-processed, individual-level cohort data used in the study. It includes variables for patient demographics, tumor characteristics, treatment history, and time-to-event outcomes.

> **Note**: This file is derived from restricted-use data and is not available in this repository.

### 2. `MortalityRatesTables_17Mar2017.csv`

This CSV file contains U.S. life table estimates used to model expected other-cause mortality in the study cohort. It is used by `OC.surv()` and other functions in `aux.fun.R` to compare disease-specific survival against population-level background mortality.

The data are derived from:

> Gangnon RE, Sprague BL, Stout NK, et al.  
> *Contribution of breast cancer to overall mortality for US women.*  
> *Medical Decision Making.* 2018;38(1_suppl):24S–31S.  
---

## License

This code is released under the [Creative Commons Attribution–NonCommercial (CC BY-NC) license](https://creativecommons.org/licenses/by-nc/4.0/).  
You are free to use and adapt the code for **non-commercial purposes** with proper attribution.

---

## Contact

For questions or collaboration inquiries, please contact:

**Marc D. Ryser**  
Duke University  
📧 [marc.ryser@duke.edu](mailto:marc.ryser@duke.edu)
