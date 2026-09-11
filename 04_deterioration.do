/*==============================================================================
  04_deterioration.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Time to glycaemic deterioration above the intensification threshold
  ($HBA1C_THRESHOLD mmol/mol).

  INPUT   cohort_long.dta
          baseline_final.dta
          censoring_dates.dta
  OUTPUT  deterioration_patient_level.dta, deterioration_analysis_ready.dta
          Cox coefficients, time-varying estimates and figures

  KEY POINTS
  ----------
  1. Event = first post-baseline HbA1c measurement at or above threshold
     (bysort patid (eventdate), so this is evaluated per patient).
  2. Exit date is capped at death, de-registration or last collection.

  NOTE ON TIME ORIGIN
  -------------------
  Follow-up starts at baseline_hba1c_date (the measured baseline), not the
  diabetes diagnosis date used in 06_mace.

  DEPENDENCY: -parmest- (ssc install parmest)
==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

log using "$logdir/04_deterioration.log", replace text

di as txt _n "{hline 78}"
di as txt "STEP 04: Time to glycaemic deterioration"
di as txt "{hline 78}"

cd "$longfmt"


/*==============================================================================
  1. BUILD THE PATIENT-LEVEL EVENT FILE

  Baseline HbA1c, BMI and their measurement dates are derived here from the
  long file
==============================================================================*/

*--- baseline HbA1c: closest measurement within +/- $BASELINE_WINDOW days ---*
use "cohort_long.dta", clear
preserve
    keep if !missing(hba1c) & !missing(eventdate) & !missing(baseline)
    keep if inrange(hba1c, $HBA1C_MIN, $HBA1C_MAX)
    gen days_from_baseline = eventdate - baseline
    keep if inrange(days_from_baseline, -$BASELINE_WINDOW, $BASELINE_WINDOW)
    gen abs_days = abs(days_from_baseline)
    bysort patid (abs_days eventdate): keep if _n == 1
    keep patid hba1c eventdate
    rename hba1c    baseline_hba1c
    rename eventdate baseline_hba1c_date
    format baseline_hba1c_date %td
    tempfile base_hba1c
    save `base_hba1c', replace
restore

*--- baseline BMI: closest measurement within +/- $BASELINE_WINDOW days ---*
preserve
    keep if !missing(bmi) & !missing(eventdate) & !missing(baseline)
    keep if inrange(bmi, $BMI_MIN, $BMI_MAX)
    gen days_from_baseline = eventdate - baseline
    keep if inrange(days_from_baseline, -$BASELINE_WINDOW, $BASELINE_WINDOW)
    gen abs_days = abs(days_from_baseline)
    bysort patid (abs_days eventdate): keep if _n == 1
    keep patid bmi
    rename bmi baseline_bmi
    tempfile base_bmi
    save `base_bmi', replace
restore

*--- keep HbA1c rows, identify threshold crossings ---*
keep if !missing(hba1c) & !missing(eventdate) & !missing(baseline)
keep if inrange(hba1c, $HBA1C_MIN, $HBA1C_MAX)

gen t_years = (eventdate - baseline) / 365.25
drop if t_years < 0
drop if t_years > $TRAJ_MAXYRS

gen byte poor = (hba1c >= $HBA1C_THRESHOLD)

/*------------------------------------------------------------------------------
  First crossing of the threshold.
------------------------------------------------------------------------------*/
bysort patid (eventdate): gen byte event = (poor == 1 & poor[_n-1] == 0) if _n > 1
bysort patid (eventdate): replace event = (poor == 1) if _n == 1

bysort patid: egen event_date = min(cond(event == 1, eventdate, .))
format event_date %td
gen byte failure = (eventdate == event_date & event == 1)

*--- attach baseline values/dates and covariates ---*
merge m:1 patid using `base_hba1c', keep(master match) nogen
merge m:1 patid using `base_bmi',   keep(master match) nogen
merge m:1 patid using "baseline_final.dta", ///
    keepusing(gender age smoking_status imd engagement death_date hes_linked ///
              htn_prev ihd_prev copd_prev highchol_prev af_prev hf_prev ckd_prev) ///
    keep(master match) nogen

drop if missing(baseline_hba1c) | missing(baseline_bmi)

quietly summarize baseline_hba1c
gen c_bhba = baseline_hba1c - r(mean)
quietly summarize baseline_bmi
gen c_bbmi = baseline_bmi - r(mean)

*--- collapse to one row per patient ---*
bysort patid: egen byte failure_patient = max(failure)
bysort patid: egen double event_date_patient = min(event_date)
format event_date_patient %td
drop failure event_date event poor
rename failure_patient   failure
rename event_date_patient event_date

bysort patid: keep if _n == 1

gen byte base_control = (baseline_hba1c < $HBA1C_THRESHOLD)

save "deterioration_patient_level.dta", replace

di as txt "Total patients: " _N
quietly count if base_control == 1
di as txt "  baseline HbA1c <$HBA1C_THRESHOLD  (analysed): " r(N)
quietly count if base_control == 0
di as txt "  baseline HbA1c >=$HBA1C_THRESHOLD (excluded): " r(N)


/*==============================================================================
  2. CENSORING AND FOLLOW-UP
==============================================================================*/

merge 1:1 patid using "$censor/censoring_dates.dta", keep(master match) nogen

capture confirm string variable regenddate
if !_rc {
    rename regenddate dereg
    gen double dereg_date = date(dereg, "DMY")
    format dereg_date %td
    drop dereg
    rename dereg_date dereg
}
capture confirm string variable lcd
if !_rc {
    gen double lcd_date = date(lcd, "DMY")
    format lcd_date %td
    drop lcd
    rename lcd_date lcd
}

*--- administrative and clinical censoring ---*
gen double study_end = baseline_hba1c_date + $FOLLOWUP_DAYS
gen double censor_date = study_end
foreach v in dereg death_date lcd {
    capture replace censor_date = min(censor_date, `v') if !missing(`v')
}
format study_end censor_date %td

*--- events after the censoring date do not count ---*
replace failure = 0 if !missing(event_date) & event_date > censor_date
gen double exit_date = cond(failure == 1, event_date, censor_date)
format exit_date %td

drop if missing(baseline_hba1c_date) | missing(exit_date)
drop if exit_date < baseline_hba1c_date

save "deterioration_analysis_ready.dta", replace


/*==============================================================================
  3. Time to first HbA1c >=$HBA1C_THRESHOLD among patients whose baseline HbA1c
  was below the threshold.
==============================================================================*/

local DET_COVARS c.age i.gender i.smoking_status ///
    i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
    i.ckd_prev i.af_prev i.hf_prev i.imd c.engagement ///
    c.c_bhba c.c_bbmi

*--- same list without c_bhba, which the tband interaction expands ---*
local DET_COVARS_NOHBA c.age i.gender i.smoking_status ///
    i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
    i.ckd_prev i.af_prev i.hf_prev i.imd c.engagement ///
    c.c_bbmi

di as txt _n "{hline 60}"
di as txt "Time to HbA1c >=$HBA1C_THRESHOLD"
di as txt "  (baseline HbA1c <$HBA1C_THRESHOLD mmol/mol)"
di as txt "{hline 60}"

use "deterioration_analysis_ready.dta", clear
keep if base_control == 1

quietly count
di as txt "  patients: " r(N)
quietly count if failure == 1
di as txt "  events:   " r(N)

stset exit_date, id(patid) origin(time baseline_hba1c_date) ///
    failure(failure == 1) scale(365.25)

*--- unadjusted KM ---*
sts graph, ///
    risktable(0(1)10, size(vsmall) title("Number at risk", size(vsmall))) ///
    xtitle("Years from baseline", margin(medium)) ///
    xlabel(0(1)10) ///
    ytitle("Patients remaining below $HBA1C_THRESHOLD mmol/mol (%)", margin(medium)) ///
    ylabel(0 "0" .2 "20" .4 "40" .6 "60" .8 "80" 1 "100", angle(horizontal)) ///
    title("") ///
    graphregion(color(white)) plotregion(color(white)) scheme(s2color)
graph export "deterioration_km_unadjusted.png", replace width(1200) height(800)

stci
sts list, at(1 2 5 9.99)

*--- adjusted Cox ---*
stcox `DET_COVARS'
estimates store cox_deterioration
estimates save "cox_deterioration.ster", replace
estat phtest, detail


/*------------------------------------------------------------------------------
  3.1 TIME-VARYING EFFECT OF BASELINE HbA1c

  The global PH test is rejected, driven almost entirely by c_bhba (rho ~ -0.21,
  contributing the large majority of the global chi2). Follow-up is split at 2
  and 5 years and c_bhba interacted with the resulting bands, giving hazard
  ratios specific to 0-2, 2-5 and 5-10 years.
------------------------------------------------------------------------------*/

preserve

    stsplit tband, at(2 5)

    di as txt _n "{hline 60}"
    di as txt "Time-varying effect of baseline HbA1c"
    di as txt "{hline 60}"

    stcox `DET_COVARS_NOHBA' i.tband##c.c_bhba
    estimates store cox_deterioration_tvc
    estimates save "cox_deterioration_tvc.ster", replace

    di as txt _n "Period-specific HR per mmol/mol higher baseline HbA1c"

    di as txt _n "  0-2 years:"
    lincom _b[c_bhba], eform

    di as txt _n "  2-5 years:"
    lincom _b[c_bhba] + _b[2.tband#c.c_bhba], eform

    di as txt _n "  5-10 years:"
    lincom _b[c_bhba] + _b[5.tband#c.c_bhba], eform

    *--- residual PH test after allowing the time-varying effect ---*
    di as txt _n "Residual PH test (split model):"
    estat phtest, detail

    capture parmest, saving("deterioration_cox_tvc_coefficients.dta", replace) norestore
    if _rc di as error "  parmest not installed - tvc coefficient file skipped"

restore


/*------------------------------------------------------------------------------
  3.2 ADJUSTED SURVIVAL CURVE
------------------------------------------------------------------------------*/

capture parmest, saving("deterioration_cox_coefficients.dta", replace) norestore
if _rc di as error "  parmest not installed - coefficient file skipped"

use "deterioration_analysis_ready.dta", clear
keep if base_control == 1

stset exit_date, id(patid) origin(time baseline_hba1c_date) ///
    failure(failure == 1) scale(365.25)
quietly stcox `DET_COVARS'

stcurve, survival range(0 10) ///
    xtitle("Years from baseline") xlabel(0(1)10) ///
    ytitle("Patients remaining below $HBA1C_THRESHOLD mmol/mol (%)", margin(medium)) ///
    ylabel(0 "0" .2 "20" .4 "40" .6 "60" .8 "80" 1 "100", angle(horizontal)) ///
    title("") ///
    graphregion(color(white)) plotregion(color(white)) scheme(s2color)
graph export "deterioration_adjusted_survival.png", replace width(1200) height(800)

log close

di as txt _n "04_deterioration complete."
