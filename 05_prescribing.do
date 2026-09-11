/*==============================================================================
  05_prescribing.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Glucose-lowering prescribing patterns.

   Section 1  Build the prescription dataset
   Section 2  Mean concurrent drug classes and proportion on each class per month
   Section 3  Treatment escalation (time to 1st, 2nd, 3rd class and gaps)
   Section 4  Time to first prescription of any class (Kaplan-Meier, full cohort)
   Section 5  Escalation by HbA1c trajectory group
   Section 6  Treatment status and timing by trajectory
   Section 7  Baseline HbA1c and engagement by prescription status
   Section 8  Predictors of ever receiving a prescription
   Section 9  Annual HbA1c change by treatment status
   Section 10 Clinical inertia (untreated despite worsening control)

  INPUT   cohort_long.dta
          baseline_demographics.dta
          hba1c_trajectories.dta
  OUTPUT  prescription_data.dta, escalation_analysis.dta,
          time_to_first_rx.dta, escalation_by_trajectory.dta


  NOTE ON DRUG2/DRUG3 DEFINITION
  ------------------------------
  drug2_date and drug3_date are the second and third DISTINCT initiation dates,
  not the second and third classes. A patient starting metformin and a
  sulfonylurea on the same day has one initiation date, so years_to_drug2 is
  the time to their next new class after that day.

==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

global RX_RECODE_DAYS 365      // pre-index prescriptions within this window are
                               // treated as initiated at diagnosis
global RX_SURVTIME    0.5      // survival time assigned to those patients (days)
global RX_COVER_DAYS  90       // assumed days of coverage per prescription issue
global RX_SENS_WINDOWS 180 90  // alternative pre-index windows for sensitivity

log using "$logdir/05_prescribing.log", replace text

di as txt _n "{hline 78}"
di as txt "STEP 05: Prescribing patterns"
di as txt "{hline 78}"

cd "$longfmt"

global GLUCOSE_CLASSES metformin su insulin glp sglt2 dppiv // TZD is not in the dataset, so it's excluded


/*==============================================================================
  1. PRESCRIPTION DATA
==============================================================================*/

use "cohort_long.dta", clear

*--- every class flag must exist; Section 2 depends on all of them ---*
local MISSING_CLASSES ""
foreach d of global GLUCOSE_CLASSES {
    capture confirm variable `d'
    if _rc local MISSING_CLASSES "`MISSING_CLASSES' `d'"
}
if "`MISSING_CLASSES'" != "" {
    di as error "Class flag(s) absent from the long file:`MISSING_CLASSES'"
    di as error "Section 2 cannot run. Check the drug code-list merge in 01/02."
    exit 111
}

gen byte has_drug = 0
foreach d of global GLUCOSE_CLASSES {
    replace has_drug = 1 if `d' == 1
}
keep if has_drug == 1

/*------------------------------------------------------------------------------
  Guard against empty results
------------------------------------------------------------------------------*/
quietly count
if r(N) == 0 {
    di as error "No prescription records retained."
    di as error "Is GLUCOSE_CLASSES populated? Current value: $GLUCOSE_CLASSES"
    di as error "prescription_data.dta has NOT been overwritten."
    exit 2000
}

local EMPTY_CLASSES ""
foreach d of global GLUCOSE_CLASSES {
    quietly count if `d' == 1
    if r(N) == 0 local EMPTY_CLASSES "`EMPTY_CLASSES' `d'"
}
if "`EMPTY_CLASSES'" != "" {
    di as error "Class flag(s) present but never set:`EMPTY_CLASSES'"
    di as error "Check the drug code-list merge in 01/02 before using these results."
    di as error "To continue without them, comment out the -exit- below."
    exit 111
}

gen double drug_date = eventdate
format drug_date %td
gen months_from_baseline = floor((drug_date - baseline) / 30.44)

save "prescription_data.dta", replace
flowcount "with >=1 glucose-lowering prescription"

/*------------------------------------------------------------------------------
  1.1 DIAGNOSTIC - do the class flags populate?
------------------------------------------------------------------------------*/

di as txt _n "Prescription records by class:"
foreach d of global GLUCOSE_CLASSES {
    quietly count if `d' == 1
    di as txt "  `d' = " r(N)
}

gen int rx_year = year(drug_date)
di as txt _n "Prescription records by class and calendar year:"
foreach d of global GLUCOSE_CLASSES {
    di as txt _n "  `d':"
    tab rx_year if `d' == 1
}
drop rx_year


/*==============================================================================
  2. Concurrent classes and class proportions by month

  Denominator is patients still under follow-up in each month, not
  prescription records.
==============================================================================*/

/*------------------------------------------------------------------------------
  2.1 RISK-SET SKELETON - one row per patient per month under follow-up
------------------------------------------------------------------------------*/

use "baseline_demographics.dta", clear
merge 1:1 patid using "$censor/censoring_dates.dta", keep(master match) nogen

capture rename regenddate dereg
foreach v in dereg lcd {
    capture confirm string variable `v'
    if !_rc {
        gen double `v'_d = date(`v', "DMY")
        format `v'_d %td
        drop `v'
        rename `v'_d `v'
    }
}

gen double fu_end = baseline + $FOLLOWUP_DAYS
foreach v in dereg death_date lcd {
    capture replace fu_end = min(fu_end, `v') if !missing(`v')
}
format fu_end %td

gen int last_month = floor((fu_end - baseline) / 30.44)
replace last_month = 120 if last_month > 120
drop if last_month < 0

keep patid baseline last_month
quietly count
di as txt _n "Patients contributing to the monthly panel: " r(N)

expand last_month + 1
bysort patid: gen int month = _n - 1

tempfile skeleton
save `skeleton', replace


/*------------------------------------------------------------------------------
  2.2 PRESCRIPTION COVERAGE - expand each issue across the months it covers
------------------------------------------------------------------------------*/

use "prescription_data.dta", clear
keep patid baseline drug_date $GLUCOSE_CLASSES

gen int m_start = floor((drug_date - baseline) / 30.44)
gen int m_end   = floor((drug_date + $RX_COVER_DAYS - baseline) / 30.44)

*--- pre-index recode, matching the policy applied in Section 4 ---*
quietly count if drug_date < baseline & drug_date >= baseline - $RX_RECODE_DAYS
di as txt "  FLOW | 4a/4b prescriptions recoded to month 0 = " r(N)
replace m_start = 0 if drug_date < baseline & drug_date >= baseline - $RX_RECODE_DAYS
replace m_end   = 0 if m_end < 0 & m_start == 0

quietly count if drug_date < baseline - $RX_RECODE_DAYS
di as txt "  FLOW | 4a/4b prescriptions dropped, >$RX_RECODE_DAYS d pre-index = " r(N)
drop if drug_date < baseline - $RX_RECODE_DAYS

replace m_end = 120 if m_end > 120
drop if m_start > 120
drop if m_start < 0

gen long rxid = _n
gen int nmonths = m_end - m_start + 1
expand nmonths
bysort rxid: gen int month = m_start + _n - 1

collapse (max) $GLUCOSE_CLASSES, by(patid month)

tempfile coverage
save `coverage', replace


/*------------------------------------------------------------------------------
  2.3 COMBINE AND SUMMARISE
------------------------------------------------------------------------------*/

use `skeleton', clear
merge 1:1 patid month using `coverage', keep(master match) nogen

foreach d of global GLUCOSE_CLASSES {
    replace `d' = 0 if missing(`d')
}

egen byte total_drugs = rowtotal($GLUCOSE_CLASSES)
gen byte any_drug = (total_drugs > 0)

*--- patients at risk per month, for the figure note and risk table ---*
preserve
    collapse (count) n_at_risk = patid, by(month)
    tempfile atrisk
    save `atrisk', replace
restore

/*------------------------------------------------------------------------------
  Conditional mean: agents per patient among those currently treated.
------------------------------------------------------------------------------*/
preserve
    collapse (mean) mean_if_treated = total_drugs if total_drugs > 0, by(month)
    tempfile condmean
    save `condmean', replace
restore

collapse (mean) total_drugs any_drug $GLUCOSE_CLASSES, by(month)
merge 1:1 month using `condmean', nogen
merge 1:1 month using `atrisk', nogen

gen years = month / 12
label variable years "Years from diagnosis"

save "prescribing_monthly_panel.dta", replace

di as txt _n "Panel summary at selected months:"
list month years n_at_risk any_drug total_drugs mean_if_treated ///
    if inlist(month, 0, 6, 12, 60, 119), noobs


/*------------------------------------------------------------------------------
  2.4 FIGURES
------------------------------------------------------------------------------*/

twoway line total_drugs years, ///
    xtitle("Years from type 2 diabetes diagnosis", margin(medium)) ///
    ytitle("Mean number of glucose-lowering agents", margin(r+6)) ///
    xlabel(0(1)10) ylabel(0(0.1)0.6, angle(horizontal)) ///
    note("Denominator includes all patients under follow-up, including those" ///
         " not currently prescribed any glucose-lowering agent.", size(tiny) position(6)) ///
    graphregion(color(white)) scheme(s2color)
graph export "concurrent_drugs_over_time.png", replace width(1200) height(800)

local plotcmd ""
foreach d of global GLUCOSE_CLASSES {
    local plotcmd "`plotcmd' (line `d' years)"
}
twoway `plotcmd', ///
    xtitle("Years from type 2 diabetes diagnosis", margin(medium)) ///
    ytitle("Proportion of patients prescribed", margin(r+6)) ///
    xlabel(0(1)10) ylabel(0(.1).5, angle(horizontal)) ///
    legend(order(1 "Metformin" 2 "SU" 3 "Insulin" 4 "GLP-1 RA" ///
                 5 "SGLT2i" 6 "DPP-4i") rows(1) position(6)) ///
    note("Denominator: patients under follow-up in each month." ///
         " A class is counted as current for $RX_COVER_DAYS days after each issue." ///
         " Thiazolidinediones are omitted: no prescriptions were recorded.", ///
         size(tiny) position(6)) ///
    graphregion(color(white)) scheme(s2color)
graph export "drug_class_proportions_over_time.png", replace width(1200) height(800)


/*==============================================================================
  3. Treatment escalation
==============================================================================*/

use "prescription_data.dta", clear

foreach d of global GLUCOSE_CLASSES {
    bysort patid: egen double `d'_first = min(cond(`d' == 1, drug_date, .))
    format `d'_first %td
}
bysort patid: keep if _n == 1

*--- explicit varlist listing each class's *_first variable ---*
local FIRSTVARS ""
foreach d of global GLUCOSE_CLASSES {
    local FIRSTVARS "`FIRSTVARS' `d'_first"
}

egen double drug1_date = rowmin(`FIRSTVARS')
gen n_classes = 0
foreach d of global GLUCOSE_CLASSES {
    replace n_classes = n_classes + (!missing(`d'_first))
}

*--- second and third DISTINCT initiation dates (see header note) ---*
gen double drug2_date = .
gen double drug3_date = .
foreach d of global GLUCOSE_CLASSES {
    replace drug2_date = min(cond(`d'_first > drug1_date, `d'_first, .), drug2_date)
}
foreach d of global GLUCOSE_CLASSES {
    replace drug3_date = min(cond(`d'_first > drug2_date, `d'_first, .), drug3_date)
}
format drug1_date drug2_date drug3_date %td

*--- how often do classes share an initiation date? ---*
quietly count if n_classes >= 2 & missing(drug2_date)
di as txt _n "  DIAG | >=2 classes but only 1 distinct initiation date = " r(N)
quietly count if n_classes >= 3 & missing(drug3_date)
di as txt "  DIAG | >=3 classes but <3 distinct initiation dates = " r(N)

*--- capture drop guards against first_drug_date already existing in the long file ---*
capture drop first_drug_date
gen double first_drug_date = drug1_date
format first_drug_date %td
gen byte ever_any_drug = !missing(first_drug_date)

gen years_to_drug1 = (drug1_date - baseline) / 365.25
gen years_to_drug2 = (drug2_date - baseline) / 365.25
gen years_to_drug3 = (drug3_date - baseline) / 365.25
gen gap_1_to_2     = years_to_drug2 - years_to_drug1
gen gap_2_to_3     = years_to_drug3 - years_to_drug2

di as txt _n "Escalation timing (years):"
tabstat years_to_drug1 years_to_drug2 years_to_drug3 gap_1_to_2 gap_2_to_3, ///
    stats(n mean sd p25 median p75) columns(statistics)

keep patid baseline first_drug_date ever_any_drug n_classes ///
     years_to_drug* gap_* `FIRSTVARS' drug?_date
save "escalation_analysis.dta", replace


/*==============================================================================
  4. Time to first prescription, whole cohort

  Never-treated patients are retained here as censored observations; the
  prescription file above contains only treated patients.
==============================================================================*/

use "baseline_demographics.dta", clear
merge 1:1 patid using "escalation_analysis.dta", ///
    keepusing(first_drug_date ever_any_drug years_to_drug1) keep(master match) nogen
replace ever_any_drug = 0 if missing(ever_any_drug)

merge 1:1 patid using "$censor/censoring_dates.dta", keep(master match) nogen
capture rename regenddate dereg
foreach v in dereg lcd {
    capture confirm string variable `v'
    if !_rc {
        gen double `v'_d = date(`v', "DMY")
        format `v'_d %td
        drop `v'
        rename `v'_d `v'
    }
}

flowcount "cohort entering Section 4 analysis"

gen byte failure = ever_any_drug

*--- censoring date: earliest of 10 years, death, de-registration, collection ---*
gen double study_end = baseline + $FOLLOWUP_DAYS
gen double censor_date = study_end
foreach v in dereg death_date lcd {
    capture replace censor_date = min(censor_date, `v') if !missing(`v')
}
format study_end censor_date %td

/*------------------------------------------------------------------------------
  Pre-index prescriptions.
  Prescription records commonly precede the clinical code for the diagnosis.
  Those within $RX_RECODE_DAYS are treated as initiated at diagnosis; those
  earlier are excluded as probable prevalent, misdated, or alternative-
  indication cases. Both counts are logged for the flow diagram.
------------------------------------------------------------------------------*/
quietly count if failure == 1 & first_drug_date <= baseline & ///
                 first_drug_date >= baseline - $RX_RECODE_DAYS
di as txt "  FLOW | prescriptions recoded to index (on index or <=$RX_RECODE_DAYS d before) = " r(N)

quietly count if failure == 1 & first_drug_date < baseline - $RX_RECODE_DAYS
di as txt "  FLOW | EXCLUDED, first prescription >$RX_RECODE_DAYS d before index = " r(N)

quietly count if censor_date < baseline
di as txt "  FLOW | EXCLUDED, censored before index = " r(N)
drop if censor_date < baseline

replace first_drug_date = baseline + $RX_SURVTIME ///
    if failure == 1 & first_drug_date <= baseline & ///
       first_drug_date >= baseline - $RX_RECODE_DAYS

drop if failure == 1 & first_drug_date < baseline

/*------------------------------------------------------------------------------
  A first prescription dated after death, de-registration or last collection
  falls outside observed follow-up. Such patients are censored rather than
  contributing an event.
------------------------------------------------------------------------------*/
quietly count if failure == 1 & first_drug_date > censor_date
di as txt "  FLOW | first prescription after censoring, recoded to censored = " r(N)
replace failure = 0 if failure == 1 & first_drug_date > censor_date

gen double exit_date = cond(failure == 1, first_drug_date, censor_date)
format exit_date %td

drop if missing(exit_date)
drop if exit_date < baseline
flowcount "Section 4 analysis cohort"

stset exit_date, id(patid) origin(time baseline) failure(failure == 1) scale(365.25)

sts graph, failure ///
    risktable(0(1)10, size(vsmall) title("Number at risk", size(vsmall))) ///
    xtitle("Years from type 2 diabetes diagnosis") ///
    ytitle("Cumulative proportion prescribed (%)") ///
    xlabel(0(1)10) ///
    ylabel(0 "0" .2 "20" .4 "40" .6 "60" .8 "80" 1 "100", angle(horizontal)) ///
    title("") ///
    note("Prescriptions within $RX_RECODE_DAYS days before diagnosis treated as initiated at diagnosis") ///
    graphregion(color(white)) plotregion(color(white)) scheme(s2color)
graph export "km_time_to_first_prescription.png", replace width(1200) height(800)

di as txt _n "Cumulative proportion prescribed:"
sts list, at(0.002 0.5 1 2 3 5 7 9.99) failure

di as txt _n "Median time to first prescription:"
stci

save "time_to_first_rx.dta", replace


/*==============================================================================
  4.1 SENSITIVITY - narrower pre-index window
==============================================================================*/

foreach w of global RX_SENS_WINDOWS {

    di as txt _n "{hline 70}"
    di as txt "SENSITIVITY: pre-index window = `w' days"
    di as txt "{hline 70}"

    preserve

        use "baseline_demographics.dta", clear
        merge 1:1 patid using "escalation_analysis.dta", ///
            keepusing(first_drug_date ever_any_drug) keep(master match) nogen
        replace ever_any_drug = 0 if missing(ever_any_drug)

        merge 1:1 patid using "$censor/censoring_dates.dta", keep(master match) nogen
        capture rename regenddate dereg
        foreach v in dereg lcd {
            capture confirm string variable `v'
            if !_rc {
                gen double `v'_d = date(`v', "DMY")
                format `v'_d %td
                drop `v'
                rename `v'_d `v'
            }
        }

        gen byte failure = ever_any_drug

        gen double censor_date = baseline + $FOLLOWUP_DAYS
        foreach v in dereg death_date lcd {
            capture replace censor_date = min(censor_date, `v') if !missing(`v')
        }
        format censor_date %td
        drop if censor_date < baseline

        quietly count if failure == 1 & first_drug_date < baseline - `w'
        di as txt "  excluded, first Rx >`w' d before index = " r(N)

        replace first_drug_date = baseline + $RX_SURVTIME ///
            if failure == 1 & first_drug_date <= baseline & ///
               first_drug_date >= baseline - `w'

        drop if failure == 1 & first_drug_date < baseline
        replace failure = 0 if failure == 1 & first_drug_date > censor_date

        gen double exit_date = cond(failure == 1, first_drug_date, censor_date)
        drop if missing(exit_date) | exit_date < baseline

        quietly count
        di as txt "  analysis cohort = " r(N)

        stset exit_date, id(patid) origin(time baseline) ///
            failure(failure == 1) scale(365.25)

        sts list, at(0.5 1 2 5 9.99) failure
        stci

    restore
}



/*==============================================================================
  5. Escalation by HbA1c trajectory group
==============================================================================*/

use "escalation_analysis.dta", clear
merge 1:1 patid using "hba1c_trajectories.dta", ///
    keepusing(trajectory b_slope indiv_slope) keep(master match) nogen
keep if !missing(trajectory)

foreach d of global GLUCOSE_CLASSES {
    gen byte ever_`d' = !missing(`d'_first)
}

di as txt _n "Time to first / second / third class, by trajectory:"
foreach v in years_to_drug1 years_to_drug2 years_to_drug3 gap_1_to_2 {
    di as txt _n "  `v':"
    tabstat `v', by(trajectory) stats(n mean sd p25 median p75)
}

di as txt _n "Proportion ever prescribed each class, by trajectory:"
tabstat ever_*, by(trajectory) stats(n mean)

save "escalation_by_trajectory.dta", replace


/*==============================================================================
  6. Treatment status and timing by trajectory
==============================================================================*/

use "time_to_first_rx.dta", clear
merge 1:1 patid using "hba1c_trajectories.dta", ///
    keepusing(trajectory b_slope indiv_slope) keep(master match) nogen

di as txt _n "Trajectory availability in the full cohort:"
tab trajectory, missing

di as txt _n "Treatment status by trajectory group:"
tab trajectory failure, row missing

stset exit_date, id(patid) origin(time baseline) failure(failure == 1) scale(365.25)

levelsof trajectory, local(TLEVELS)
local TLEGEND ""
local i = 0
foreach t of local TLEVELS {
    local ++i
    local tlab : label (trajectory) `t'
    local TLEGEND `"`TLEGEND' `i' "`tlab'""'
}

sts graph, failure by(trajectory) ///
    risktable(0(2)10, size(vsmall) title("Number at risk", size(vsmall))) ///
    xtitle("Years from type 2 diabetes diagnosis") ///
    ytitle("Cumulative proportion prescribed (%)") ///
    xlabel(0(1)10) ///
    ylabel(0 "0" .2 "20" .4 "40" .6 "60" .8 "80" 1 "100", angle(horizontal)) ///
    legend(order(`TLEGEND') rows(1) position(6)) ///
    title("") ///
    graphregion(color(white)) plotregion(color(white)) scheme(s2color)
graph export "km_prescribing_by_trajectory.png", replace width(1200) height(800)

sts test trajectory
sts list, at(1 2 5 9.99) by(trajectory) failure

*--- untreated at 10 years among worsening patients above threshold ---*
capture confirm variable baseline_hba1c
if _rc {
    merge 1:1 patid using "baseline_final.dta", ///
        keepusing(baseline_hba1c) keep(master match) nogen
}

*--- confirm which level is the worsening group before selecting on it ---*
local WORSENING 3
local wlab : label (trajectory) `WORSENING'
di as txt _n "Selecting trajectory == `WORSENING', labelled: `wlab'"
di as txt "Check this is the worsening group before reporting the counts below."

di as txt _n "Untreated within 10y, worsening AND baseline HbA1c >=$HBA1C_THRESHOLD:"
count if trajectory == `WORSENING' & baseline_hba1c >= $HBA1C_THRESHOLD
count if trajectory == `WORSENING' & baseline_hba1c >= $HBA1C_THRESHOLD & failure == 0

save "treatment_by_trajectory.dta", replace



/*==============================================================================
  7. BASELINE HbA1c AND ENGAGEMENT BY PRESCRIPTION STATUS

  Restricted to patients with a recorded baseline HbA1c.
==============================================================================*/

di as txt _n "{hline 70}"
di as txt "Baseline HbA1c and engagement by prescription status"
di as txt "{hline 70}"

use "time_to_first_rx.dta", clear
keep patid baseline ever_any_drug failure engagement

/*------------------------------------------------------------------------------
  baseline_hba1c is derived in 04_deterioration.do from the long file, not
  merged from baseline_final
------------------------------------------------------------------------------*/
merge 1:1 patid using "deterioration_patient_level.dta", ///
    keepusing(baseline_hba1c baseline_bmi) keep(master match) nogen

quietly count
di as txt _n "Section 4 cohort: " r(N)
quietly count if !missing(baseline_hba1c)
di as txt "With recorded baseline HbA1c: " r(N)

keep if !missing(baseline_hba1c)

gen byte ever_rx = ever_any_drug
label define rxlab 0 "Never prescribed" 1 "Ever prescribed"
label values ever_rx rxlab
tab ever_rx, missing

*--- HbA1c categories as reported ---*
gen byte hba1c_cat = .
replace hba1c_cat = 1 if baseline_hba1c < 48
replace hba1c_cat = 2 if inrange(baseline_hba1c, 48, 52)
replace hba1c_cat = 3 if inrange(baseline_hba1c, 53, 57)
replace hba1c_cat = 4 if inrange(baseline_hba1c, 58, 74)
replace hba1c_cat = 5 if baseline_hba1c >= 75
label define hbacat 1 "<48" 2 "48-52" 3 "53-57" 4 "58-74" 5 ">=75"
label values hba1c_cat hbacat

gen byte below_53 = (baseline_hba1c < 53)

di as txt _n "Baseline HbA1c, mean (SD) by prescription status:"
tabstat baseline_hba1c, by(ever_rx) stats(n mean sd) format(%9.1f)
ttest baseline_hba1c, by(ever_rx)
regress baseline_hba1c i.ever_rx

di as txt _n "Healthcare engagement, mean contacts (SD) by prescription status:"
tabstat engagement, by(ever_rx) stats(n mean sd) format(%9.1f)
ttest engagement, by(ever_rx)

di as txt _n "HbA1c category by prescription status:"
tab hba1c_cat ever_rx, col chi2

di as txt _n "Below 53 mmol/mol by prescription status:"
tab below_53 ever_rx, col chi2

*--- sub-diagnostic-threshold values, both groups ---*
di as txt _n "Baseline HbA1c <48 mmol/mol, by prescription status:"
forvalues g = 0/1 {
    quietly count if ever_rx == `g'
    local den = r(N)
    quietly count if ever_rx == `g' & baseline_hba1c < 48
    local lab : label rxlab `g'
    di as txt "  " %-18s "`lab'" r(N) " of `den' (" %4.1f 100*r(N)/`den' "%)"
}

save "prescription_status_by_hba1c.dta", replace


/*==============================================================================
  8. PREDICTORS OF EVER RECEIVING A PRESCRIPTION
==============================================================================*/

di as txt _n "{hline 70}"
di as txt "Predictors of ever receiving a prescription"
di as txt "{hline 70}"

use "time_to_first_rx.dta", clear
keep patid baseline ever_any_drug engagement

merge 1:1 patid using "deterioration_patient_level.dta", ///
    keepusing(baseline_hba1c baseline_bmi) keep(master match) nogen

merge 1:1 patid using "baseline_final.dta", ///
    keepusing(gender age smoking_status imd htn_prev ihd_prev copd_prev ///
              highchol_prev ckd_prev af_prev hf_prev) ///
    keep(master match) nogen

local RX_COVARS c.age i.gender i.smoking_status ///
    i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
    i.ckd_prev i.af_prev i.hf_prev i.imd c.engagement ///
    c.baseline_hba1c c.baseline_bmi

logistic ever_any_drug `RX_COVARS'
estimates store logit_m5
estimates save "logit_ever_prescribed.ster", replace

di as txt _n "Complete-case n: " e(N)

/*------------------------------------------------------------------------------
  per-10 mmol/mol scaling
------------------------------------------------------------------------------*/
di as txt _n "Baseline HbA1c, per 10 mmol/mol:"
lincom 10*_b[baseline_hba1c], or

di as txt _n "Age, per decade:"
lincom 10*_b[age], or

di as txt _n "Discrimination:"
lroc, nograph
capture estat gof, group(10)
if _rc di as error "  estat gof unavailable"

capture parmest, eform saving("logistic_ever_prescribed_coefficients.dta", replace) norestore
if _rc di as error "  parmest not installed - coefficient file skipped"


/*==============================================================================
  9. ANNUAL HbA1c CHANGE, TREATED VERSUS UNTREATED

  indiv_slope is the per-patient annual change in HbA1c estimated in
  03_trajectories.do.
==============================================================================*/

di as txt _n "{hline 70}"
di as txt "ANNUAL HbA1c CHANGE BY TREATMENT STATUS"
di as txt "{hline 70}"

use "time_to_first_rx.dta", clear
keep patid ever_any_drug

merge 1:1 patid using "hba1c_trajectories.dta", ///
    keepusing(trajectory indiv_slope) keep(master match) nogen

keep if !missing(indiv_slope)

quietly count
di as txt _n "Patients with an estimable slope: " r(N)

label values ever_any_drug rxlab

di as txt _n "Annual HbA1c change (mmol/mol/year) by treatment status:"
tabstat indiv_slope, by(ever_any_drug) stats(n mean sd p25 median p75) format(%9.2f)
ttest indiv_slope, by(ever_any_drug)
regress indiv_slope i.ever_any_drug

di as txt _n "Trajectory group by treatment status:"
tab trajectory ever_any_drug, col chi2


/*==============================================================================
 10. CLINICAL INERTIA

  Patients with a worsening trajectory who were already at or above the
  intensification threshold at baseline. Both the crude proportion and the
  censoring-adjusted estimate
==============================================================================*/

di as txt _n "{hline 70}"
di as txt "CLINICAL INERTIA: worsening trajectory, baseline HbA1c >=$HBA1C_THRESHOLD"
di as txt "{hline 70}"

use "treatment_by_trajectory.dta", clear

capture confirm variable baseline_hba1c
if _rc {
    merge 1:1 patid using "deterioration_patient_level.dta", ///
        keepusing(baseline_hba1c) keep(master match) nogen
}

local WORSENING 3
local wlab : label (trajectory) `WORSENING'
di as txt _n "trajectory == `WORSENING' is labelled: `wlab'"

quietly count if trajectory == `WORSENING' & baseline_hba1c >= $HBA1C_THRESHOLD
local n_denom = r(N)
quietly count if trajectory == `WORSENING' & baseline_hba1c >= $HBA1C_THRESHOLD & failure == 0
local n_untreated = r(N)

di as txt _n "Crude: `n_untreated' of `n_denom' (" %4.1f 100*`n_untreated'/`n_denom' "%) untreated"

preserve
    keep if trajectory == `WORSENING' & baseline_hba1c >= $HBA1C_THRESHOLD
    stset exit_date, id(patid) origin(time baseline) failure(failure == 1) scale(365.25)

    di as txt _n "Censoring-adjusted cumulative proportion prescribed:"
    sts list, at(1 2 5 9.99) failure

    /*--------------------------------------------------------------------------
      Follow-up among the untreated. This is observed time WITHOUT treatment,
      not time to prescription - these patients never received one. It tests
      whether the untreated proportion reflects sustained non-treatment or
      simply short observation.
    --------------------------------------------------------------------------*/
    di as txt _n "Observed follow-up among the untreated (years):"
    summarize _t if failure == 0, detail
restore


log close

di as txt _n "05_prescribing complete."
