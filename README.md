# EPI-DIABESITY
Stata code for a study of longitudinal clinical outcomes among adults living with obesity at the time of type 2 diabetes diagnosis, using CPRD Aurum primary care data

## Data availability

CPRD data cannot be shared or re-hosted, so no data files are included in this repository. `00_config.do` uses placeholder paths throughout. Edit these to point at your own local copy of the extract before running anything else.

## Requirements

- Stata 18 or later
- `parmest` (`ssc install parmest`), used in `04_deterioration.do` and `06_mace.do`

## Pipeline

Run the scripts in order. Each one begins with `if "$root" == "" do "00_config.do"`, so `00_config.do` only needs to be run once per session. Every other script can also be run standalone, since it will load the config itself if `$root` isn't already set.

| Script | Builds | Depends on |
|---|---|---|
| `00_config.do` | Shared paths, cohort/cleaning thresholds, code lists, and the covariate adjustment set. Loads no data itself. | — |
| `01_cohort_cleaning.do` | Cohort eligibility. | 00 |
| `02_baseline_demographics.do` | One-row-per-patient baseline file and baseline summary statistics. Adds ethnicity/IMD, derives baseline HbA1c/BMI, splits comorbidities into prevalent/incident, and derives `hes_linked`, the HES/ONS linkage-eligibility flag used to restrict the MACE analyses in 06. | 00, 01 |
| `03_trajectories.do` | Longitudinal HbA1c and BMI trajectories: population-averaged marginal means, individual-level trajectory estimates (improving/stable/worsening, relative to the cohort-average slope), and annual rate-of-change by subgroup (age band, sex, IMD tertile, ethnicity). | 00, 01, 02 |
| `04_deterioration.do` | Time to glycaemic deterioration above the intensification threshold, among patients with baseline HbA1c below it. Cox model plus a time-varying-effect variant for baseline HbA1c. | 00, 01, 02 |
| `05_prescribing.do` | Glucose-lowering prescribing patterns: concurrent drug classes and class proportions over time, treatment escalation (time to 1st/2nd/3rd distinct drug class), time to first prescription (Kaplan-Meier, overall and by trajectory group), and predictors of ever being prescribed. | 00, 01, 02, 03 |
| `06_mace.do` | MACE-5 composite outcome (ACS, stroke, heart failure, TIA, CV death) built from HES + ONS, and all associated analyses: primary Cox model, MACE by trajectory group, competing-risk cumulative incidence, missing-data characterisation, multiple imputation, trajectory-threshold sensitivity, and a time-varying prescribing exposure. | 00, 01, 02, 03, 05 |
