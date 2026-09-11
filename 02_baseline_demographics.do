/*==============================================================================
  02_baseline_demographics.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Builds the one-row-per-patient baseline file and produces baseline summary
  statistics.

  INPUT   cohort_long.dta   
          ethnicity_lookup.dta, imd_lookup.dta 
  OUTPUT  cohort_long.dta   
          baseline_final.dta
          $logdir/02_baseline_summary.log
==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

di as txt _n "{hline 78}"
di as txt "STEP 02: Baseline demographics"
di as txt "{hline 78}"

cd "$longfmt"


/*==============================================================================
  1. MERGE ETHNICITY AND IMD INTO THE LONG FILE

  ethnicity_lookup.dta and imd_lookup.dta are prepared from CPRD's HES-linked
  patient reference data (the hes_patient and patient_2019_imd linkage files).
==============================================================================*/

use "cohort_long.dta", clear

rename earliest_t2d baseline
format baseline %td

* age recomputed from the single index date so it is constant within patient
capture drop age_at_t2d
gen age = year(baseline) - yob

merge m:1 patid using "$refdata/ethnicity_lookup.dta", ///
    keep(master match) nogen
merge m:1 patid using "$refdata/imd_lookup.dta", ///
    keep(master match) nogen
rename e2019_imd_10 imd

/*------------------------------------------------------------------------------
  HES/ONS linkage eligibility.
  IMD is supplied only for linkage-eligible patients
------------------------------------------------------------------------------*/
gen byte hes_linked = !missing(imd)
label define linklab 0 "Not linkage-eligible" 1 "Linkage-eligible", replace
label values hes_linked linklab

save "cohort_long.dta", replace
flowcount "long file with ethnicity and IMD"


/*==============================================================================
  2. DEMOGRAPHICS  (one row per patient)
==============================================================================*/

bysort patid: keep if _n == 1
keep patid baseline age gender yob ethnicity imd hes_linked ///
     smoking_status engagement death death_date
save "baseline_demographics.dta", replace
di as txt "Demographics: " _N " patients"


/*==============================================================================
  3. BASELINE BMI AND HbA1c
     Closest measurement to the index date within +/- $BASELINE_WINDOW days;
     ties broken by taking the earlier measurement.
==============================================================================*/

foreach m in bmi hba1c {

    use "cohort_long.dta", clear
    keep if !missing(`m') & !missing(eventdate) & !missing(baseline)

    if "`m'" == "bmi"   keep if inrange(bmi,   $BMI_MIN,   $BMI_MAX)
    if "`m'" == "hba1c" keep if inrange(hba1c, $HBA1C_MIN, $HBA1C_MAX)

    gen days_from_baseline = eventdate - baseline
    keep if inrange(days_from_baseline, -$BASELINE_WINDOW, $BASELINE_WINDOW)

    gen abs_days = abs(days_from_baseline)
    bysort patid (abs_days eventdate): keep if _n == 1

    keep patid `m' eventdate
    rename eventdate `m'_date
    rename `m' baseline_`m'

    save "baseline_`m'.dta", replace
    di as txt "Baseline `m': " _N " patients"
}


/*==============================================================================
  4. COMORBIDITIES  (prevalent at/before index vs incident after)
==============================================================================*/

use "cohort_long.dta", clear

local conds htn ihd hf pvd copd highchol ckd af depression nafld

foreach v of local conds {
    capture confirm variable `v'
    if _rc {
        di as txt "`v' not found - skipping"
        continue
    }
    bysort patid: egen `v'_prev_date = min(cond(`v'==1 & eventdate <  baseline, eventdate, .))
    bysort patid: egen `v'_post_date = min(cond(`v'==1 & eventdate >= baseline, eventdate, .))

    gen byte `v'_prev = !missing(`v'_prev_date)
    gen byte `v'_inc  = missing(`v'_prev_date) & !missing(`v'_post_date)
    gen `v'_date_final = cond(!missing(`v'_prev_date), `v'_prev_date, `v'_post_date)
    format `v'_prev_date `v'_post_date `v'_date_final %td
}

bysort patid: keep if _n == 1
keep patid *_prev *_inc *_date_final
save "baseline_comorbidities.dta", replace


/*==============================================================================
  5. ASSEMBLE THE BASELINE FILE
==============================================================================*/

use "baseline_demographics.dta", clear
foreach f in bmi hba1c comorbidities {
    merge 1:1 patid using "baseline_`f'.dta", keep(master match) nogen
}

compress
save "baseline_final.dta", replace

di as txt _n "Baseline file assembled."
count
di as txt "  linkage-eligible: "
count if hes_linked == 1
di as txt "  with baseline HbA1c: "
count if !missing(baseline_hba1c)
di as txt "  with baseline BMI: "
count if !missing(baseline_bmi)


/*==============================================================================
  6. BASELINE SUMMARY STATISTICS
==============================================================================*/

use "baseline_final.dta", clear
quietly count
local total_n = r(N)

capture log close
log using "$logdir/02_baseline_summary.log", replace text

di _n "==================================================="
di "BASELINE CHARACTERISTICS (N = `total_n')"
di "==================================================="

*--- continuous variables: mean (SD), median (IQR), missing ---*
di _n "----------------------------------------------------"
di "CONTINUOUS VARIABLES"
di "----------------------------------------------------"

foreach v in age baseline_hba1c baseline_bmi engagement {
    capture confirm variable `v'
    if _rc {
        di %-34s "`v'" "not found"
        continue
    }
    quietly summarize `v', detail
    local m  = string(r(mean), "%4.1f")
    local sd = string(r(sd),   "%4.1f")
    local p50 = string(r(p50), "%4.1f")
    local p25 = string(r(p25), "%4.1f")
    local p75 = string(r(p75), "%4.1f")
    quietly count if missing(`v')
    local miss = r(N)
    quietly count if !missing(`v')
    local n = r(N)
    di _n "`v':"
    di "  N with measurement:            " `n'
    di "  Mean (SD):                     `m' (`sd')"
    di "  Median (IQR):                  `p50' (`p25' - `p75')"
    di "  Missing:                       " `miss'
}

*--- categorical variables: level counts and percentages ---*
di _n "----------------------------------------------------"
di "CATEGORICAL VARIABLES"
di "----------------------------------------------------"

foreach v in gender ethnicity imd smoking_status hes_linked {
    capture confirm variable `v'
    if _rc {
        di %-34s "`v'" "not found"
        continue
    }
    di _n "`v':"
    quietly levelsof `v', local(lv)
    foreach l of local lv {
        quietly count if `v' == `l'
        local n = r(N)
        local pct = string(100 * `n' / `total_n', "%4.1f")
        di "  Level `l', n (%):               " `n' " (`pct'%)"
    }
    quietly count if missing(`v')
    di "  Missing:                       " r(N)
}

*--- prevalent comorbidities ---*
di _n "----------------------------------------------------"
di "CO-MORBIDITIES (prevalent at or before baseline)"
di "----------------------------------------------------"

foreach v in htn ihd hf copd highchol ckd af pvd depression nafld {
    capture confirm variable `v'_prev
    if _rc {
        di %-34s "`v'" "not found"
        continue
    }
    quietly count if `v'_prev == 1
    local n = r(N)
    local pct = string(100 * `n' / `total_n', "%4.1f")
    di %-34s "`v', n (%):" " " `n' " (`pct'%)"
}

di _n "----------------------------------------------------"
di "SUMMARY"
di "----------------------------------------------------"
di "Total patients:                    `total_n'"
foreach v in baseline_bmi baseline_hba1c {
    quietly count if !missing(`v')
    di %-34s "With `v':" " " r(N)
}
quietly count if hes_linked == 1
di %-34s "Linkage-eligible:" " " r(N)
di "----------------------------------------------------"

log close

di as txt _n "02_baseline_demographics complete."
di as txt "Baseline summary statistics written to $logdir/02_baseline_summary.log"
