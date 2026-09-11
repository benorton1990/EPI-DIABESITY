/*==============================================================================
  00_config.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Shared configuration. Every other script begins with:

      if "$root" == "" do "00_config.do"

  This file defines the paths, adjustment set, code lists, and cleaning
  thresholds.
==============================================================================*/

clear all
set more off
set linesize 100
version 18


/*------------------------------------------------------------------------------
  PATHS
  EDIT THESE to point at your own local copy of the CPRD Aurum extract
------------------------------------------------------------------------------*/
global root     "<path to your CPRD Aurum data extract>"
global obsdir   "$root/raw_extract_part_1/observation"
global longfmt  "$root/processed/cleaned_long_format"
global linked   "$root/processed/linked/final"
global refdata  "$root/processed/reference_data"      // ethnicity, IMD lookups
global censor   "$root"
global consdir  "<path to your CPRD consultations extract>"
global logdir   "$root/logs"

capture mkdir "$logdir"


/*------------------------------------------------------------------------------
  COHORT DEFINITION
------------------------------------------------------------------------------*/
global STUDY_START   "01/01/2005"   // earliest eligible T2D diagnosis
global STUDY_END     "31/12/2023"   // latest eligible T2D diagnosis
global FOLLOWUP_DAYS 3652           // 10 years
global MIN_AGE       18

* Exclusion windows (days relative to index date; negative = before)
global BARIATRIC_WINDOW  90     // bariatric surgery >90d before index excluded
global PREGNANCY_WINDOW  365    // pregnancy within +/-365d excluded
global OBESITY_WINDOW    365    // obesity code >365d after index excluded
global CANCER_WINDOW     730    // solid cancer within +/-730d excluded
global MIN_REGISTRATION  42     // >=6 weeks registered before index


/*------------------------------------------------------------------------------
  MEASUREMENT CLEANING  (applied identically wherever these appear)
------------------------------------------------------------------------------*/
global HBA1C_MIN   20
global HBA1C_MAX   200
global BMI_MIN     10
global BMI_MAX     80
global BASELINE_WINDOW 90       // +/- days around index for baseline value
global HBA1C_THRESHOLD 58       // intensification threshold (mmol/mol)


/*------------------------------------------------------------------------------
  TRAJECTORY CLASSIFICATION
------------------------------------------------------------------------------*/
global TRAJ_CUT    1            // +/- 1 unit/year around the cohort average
global TRAJ_MAXYRS 10           // longitudinal observations restricted to 0-10y


/*------------------------------------------------------------------------------
  MACE-5 DEFINITION
------------------------------------------------------------------------------*/
global ACS_PREFIX   "I21 I22"
global ACS_EXACT    "I20.0"
global STROKE_CODES "I61 I63 I64"
global HF_CODES     "I50"
global TIA_CODES    "G45"
global TIA_EXCLUDE_ICD "G45.4"
global CVDEATH_VAR  "cv_death"                // underlying cause, ICD-10 I00-I99
global PREVALENT_MACE_DAYS 90                 // events <90d after index = prevalent
global CVDEATH_EXCLUDE_ICD10 "I80 I81 I82 I83 I84 I85 I86 I87 I88 I89 I95 I99" // strict definition of cv_death

/*------------------------------------------------------------------------------
  ADJUSTMENT SET
------------------------------------------------------------------------------*/
global COVARS i.gender c.age i.smoking_status ///
    i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
    i.ckd_prev i.af_prev i.hf_prev i.imd


/*------------------------------------------------------------------------------
  MULTIPLE IMPUTATION
------------------------------------------------------------------------------*/
global MI_N      40
global MI_SEED   1234
global MI_BURNIN 10


/*------------------------------------------------------------------------------
  HELPER: log a patient count for the study flow diagram
  Usage:  flowcount "after age restriction"
------------------------------------------------------------------------------*/
capture program drop flowcount
program define flowcount
    args label
    quietly {
        tempvar tag
        egen `tag' = tag(patid)
        count if `tag'
        local n = r(N)
        drop `tag'
    }
    noisily di as txt "  FLOW | " %-45s "`label'" " | patients = " as res `n'
end


di as txt _n "{hline 78}"
di as txt "Config loaded."
di as txt "  Study window : $STUDY_START to $STUDY_END, $FOLLOWUP_DAYS days follow-up"
di as txt "  HbA1c range  : $HBA1C_MIN-$HBA1C_MAX mmol/mol"
di as txt "  BMI range    : $BMI_MIN-$BMI_MAX kg/m2"
di as txt "  MACE-5       : ACS($ACS_PREFIX $ACS_EXACT) Stroke($STROKE_CODES)"
di as txt "                 HF($HF_CODES) TIA($TIA_CODES) CVdeath($CVDEATH_VAR)"
di as txt "{hline 78}" _n
