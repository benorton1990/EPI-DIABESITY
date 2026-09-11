/*==============================================================================
  06_mace.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  MACE-5 outcomes and all associated sensitivity analyses.

    Section 1   Build the MACE-5 event dataset (HES principal diagnosis + ONS CV death)
    Section 2   Build the analytic spine, apply the linkage restriction
    Section 3   Flag prevalent MACE
    Section 4   Unadjusted incidence and Kaplan-Meier
    Section 5   Primary complete-case Cox
    Section 6   MACE by HbA1c / BMI / combined trajectory, and continuous slopes
    Section 7   Competing-risk cumulative incidence
    Section 8   Missing data characterisation
    Section 9   Multiple imputation
    Section 10  Trajectory threshold sensitivity
    Section 11  Time-varying prescribing sensitivity

  INPUT   raw_hes_diagnosis_extract.dta, raw_ons_death_registration.dta  (raw)
          baseline_final.dta               
          hba1c_trajectories.dta, bmi_trajectories.dta   
          escalation_analysis.dta                          
          censoring_dates.dta
  OUTPUT  mace_events.dta, mace_analytic_full.dta,
          mace_incident_analysis.dta, mace_mi_imputed.dta

==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

log using "$logdir/06_mace.log", replace text

di as txt _n "{hline 78}"
di as txt "STEP 06: MACE-5 outcomes"
di as txt "  ACS       : $ACS_PREFIX  $ACS_EXACT"
di as txt "  Stroke    : $STROKE_CODES"
di as txt "  HF        : $HF_CODES"
di as txt "  TIA       : $TIA_CODES  (excl. $TIA_EXCLUDE_ICD)"
di as txt "  CV death  : I-chapter, excl. $CVDEATH_EXCLUDE_ICD10"
di as txt "{hline 78}"


/*==============================================================================
  1. BUILD THE MACE-5 EVENT DATASET
==============================================================================*/

cd "$linked"

*--- 1.1: HES hospitalisations, principal diagnosis only ---*
use "raw_hes_diagnosis_extract.dta", clear
keep if d_order == 1

gen byte acs_mi = 0
foreach c of global ACS_PREFIX {
    replace acs_mi = 1 if substr(icd,1,3) == "`c'"
}
gen byte acs_ua = 0
foreach c of global ACS_EXACT {
    replace acs_ua = 1 if icd == "`c'"
}
gen byte acs = (acs_mi == 1 | acs_ua == 1)

gen byte stroke = 0
foreach c of global STROKE_CODES {
    replace stroke = 1 if substr(icd,1,3) == "`c'"
}

gen byte hf = 0
foreach c of global HF_CODES {
    replace hf = 1 if substr(icd,1,3) == "`c'"
}

gen byte tia = 0
foreach c of global TIA_CODES {
    replace tia = 1 if substr(icd,1,3) == "`c'"
}
foreach c of global TIA_EXCLUDE_ICD {
    replace tia = 0 if icd == "`c'"
}

*--- fail loudly rather than silently producing an empty component ---*
di as txt _n "Principal-diagnosis record counts:"
foreach pair in "acs_mi ACS_PREFIX" "acs_ua ACS_EXACT" "stroke STROKE_CODES" ///
                "hf HF_CODES" "tia TIA_CODES" {
    local v : word 1 of `pair'
    local g : word 2 of `pair'
    quietly count if `v' == 1
    di as txt "  `v': " r(N)
    if r(N) == 0 & trim("${`g'}") != "" {
        di as error "STOP: `v' matched zero records but codes were specified."
        di as error "HES codes carry decimals (I20.0) - check the match type."
        exit 459
    }
}
keep if acs == 1 | stroke == 1 | hf == 1 | tia == 1

gen double epistart_date = daily(epistart, "YMD")
format epistart_date %td

*--- component priority and ACS sub-priority, for same-date ties ---*
gen byte prio  = cond(acs, 1, cond(stroke, 2, cond(hf, 3, 4)))
gen byte prio2 = cond(acs_mi, 1, 2)

gen str15 mace_type = ""
replace mace_type = "ACS"           if prio == 1
replace mace_type = "Stroke"        if prio == 2
replace mace_type = "Heart failure" if prio == 3
replace mace_type = "TIA"           if prio == 4

gen str15 acs_subtype = ""
replace acs_subtype = "MI"              if prio == 1 & acs_mi == 1
replace acs_subtype = "Unstable angina" if prio == 1 & acs_mi == 0 & acs_ua == 1

rename epistart_date mace_hosp_date
bysort patid (mace_hosp_date prio prio2): keep if _n == 1

keep patid mace_hosp_date mace_type acs_subtype acs acs_mi acs_ua stroke hf tia
rename mace_hosp_date first_hosp_mace_date
tempfile hes_mace
save `hes_mace', replace
di as txt "HES MACE patients: " _N

*--- 1.2: ONS cardiovascular death, built from the raw death extract.
use "raw_ons_death_registration.dta", clear

capture confirm variable s_underlying_cod_icd10
if _rc {
    di as error "STOP: s_underlying_cod_icd10 not found — check variable name in this extract."
    exit 111
}
capture confirm variable reg_date_of_death
if _rc {
    di as error "STOP: reg_date_of_death not found — check variable name in this extract."
    exit 111
}

*--- death date: raw field is str10, "YYYY-MM-DD" 
gen double death_date = daily(reg_date_of_death, "YMD")
format death_date %td

quietly count if !missing(reg_date_of_death) & missing(death_date)
if r(N) > 0 {
    di as error "WARNING: " r(N) " non-missing reg_date_of_death failed to parse as YMD — check for stray formats."
}

*--- cv_death: any I-chapter underlying cause, excluding venous/lymphatic
*    (I80-I89), hypotension (I95), and unspecified (I99) ---*
gen byte cv_death = 0
replace cv_death = 1 if substr(s_underlying_cod_icd10, 1, 1) == "I"

quietly count if cv_death == 1
di as txt _n "CV deaths, before exclusions: " r(N)
if r(N) == 0 {
    di as error "STOP: cv_death matched zero records — check ICD-10 field format."
    exit 459
}

foreach c of global CVDEATH_EXCLUDE_ICD10 {
    replace cv_death = 0 if substr(s_underlying_cod_icd10, 1, 3) == "`c'"
}

quietly count if cv_death == 1
di as txt "CV deaths, final (excl. $CVDEATH_EXCLUDE_ICD10): " r(N)
di as txt "  Breakdown of excluded codes:"
tab s_underlying_cod_icd10 if substr(s_underlying_cod_icd10,1,1) == "I" & cv_death == 0, sort

gen double cv_death_date = death_date if cv_death == 1
format cv_death_date %td

collapse (min) cv_death_date death_date, by(patid)
format cv_death_date death_date %td

quietly count if !missing(cv_death_date)
di as txt _n "Patients with a CV death, after collapse: " r(N)

tempfile cv_death
save `cv_death', replace

*--- 1.3: combine ---*
use `hes_mace', clear
merge 1:1 patid using `cv_death', nogen

egen double mace5_date = rowmin(first_hosp_mace_date cv_death_date)
format mace5_date %td

gen byte cv_death_first = 0
replace cv_death_first = 1 if !missing(cv_death_date) & missing(first_hosp_mace_date)
replace cv_death_first = 1 if !missing(cv_death_date) & !missing(first_hosp_mace_date) ///
                            & cv_death_date <= first_hosp_mace_date

gen byte mace5_event = !missing(mace5_date)
replace mace_type   = "CV death" if cv_death_first == 1
replace acs_subtype = ""         if cv_death_first == 1

foreach v in acs acs_mi acs_ua stroke hf tia {
    replace `v' = 0 if missing(`v')
}

save "$linked/mace_events.dta", replace
di as txt _n "MACE-5 index event types:"
tab mace_type if mace5_event == 1


/*==============================================================================
  2. ANALYTIC SPINE
==============================================================================*/

cd "$longfmt"
use "baseline_final.dta", clear
flowcount "baseline cohort"

*--- 2.1: HES/ONS linkage restriction (cardiovascular analyses only) ---*
keep if hes_linked == 1
flowcount "linkage-eligible (MACE analysis cohort)"

*--- 2.2: prescribing ---*
merge 1:1 patid using "escalation_analysis.dta", ///
    keepusing(years_to_drug1 first_drug_date ever_any_drug) keep(master match) nogen
replace ever_any_drug = 0 if missing(ever_any_drug)

*--- 2.3: trajectories ---*
merge 1:1 patid using "hba1c_trajectories.dta", ///
    keepusing(trajectory b_slope indiv_slope) keep(master match) nogen
rename (trajectory b_slope indiv_slope) (hba1c_trajectory hba1c_b_slope hba1c_indiv_slope)

merge 1:1 patid using "bmi_trajectories.dta", ///
    keepusing(trajectory b_slope indiv_slope) keep(master match) nogen
rename (trajectory b_slope indiv_slope) (bmi_trajectory bmi_b_slope bmi_indiv_slope)

*--- 2.4: MACE events ---*
merge 1:1 patid using "$linked/mace_events.dta", ///
    keepusing(mace5_date mace5_event mace_type acs_subtype acs acs_mi acs_ua ///
              stroke hf tia cv_death_date cv_death_first death_date) ///
    keep(master match) nogen
replace mace5_event = 0 if missing(mace5_event)

* components must be 0, not missing, or component -stset- calls silently
* exclude observations and each rate gets a different denominator
foreach v in acs acs_mi acs_ua stroke hf tia cv_death_first {
    replace `v' = 0 if missing(`v')
}

*--- 2.5: censoring ---*
merge 1:1 patid using "$censor/censoring_dates.dta", keep(master match) nogen
capture rename regenddate dereg
confirm variable dereg          // errors if the rename did not fire

foreach v in dereg lcd mace5_date death_date {
    capture confirm string variable `v'
    if !_rc {
        gen double `v'_d = date(`v', "DMY")
        format `v'_d %td
        drop `v'
        rename `v'_d `v'
    }
}

*--- 2.6: time variables ---*
gen double diag_date = baseline
format diag_date %td
gen diagnosis_year = year(diag_date)
gen days_to_mace   = mace5_date - diag_date if !missing(mace5_date)
gen years_to_mace  = days_to_mace / 365.25

*--- 2.7: exit date and failure ---*
gen double study_end = diag_date + $FOLLOWUP_DAYS
gen double exit_date = study_end
foreach v in death_date dereg lcd mace5_date {
    replace exit_date = min(exit_date, `v') if !missing(`v')
}
format exit_date study_end %td

gen byte failure = 0
replace failure = 1 if !missing(mace5_date) & mace5_date >= diag_date & mace5_date <= exit_date

quietly count if !missing(mace5_date) & mace5_date < diag_date
di as txt "  FLOW | EXCLUDED, MACE before index date = " r(N)

drop if missing(diag_date) | missing(exit_date)
drop if exit_date < diag_date
flowcount "after removing pre-index MACE and invalid dates"

*--- 2.8: centred baseline covariates ---*
quietly summarize baseline_hba1c
gen c_hba = baseline_hba1c - r(mean)
quietly summarize baseline_bmi
gen c_bmi = baseline_bmi - r(mean)

save "mace_analytic_full.dta", replace


/*==============================================================================
  3. PREVALENT MACE
==============================================================================*/

cd "$longfmt"
use "mace_analytic_full.dta", clear
gen byte prevalent_mace = (days_to_mace < $PREVALENT_MACE_DAYS) if !missing(days_to_mace)
replace prevalent_mace = 0 if missing(prevalent_mace)
quietly count if prevalent_mace == 1
di as txt "  FLOW | EXCLUDED, MACE within $PREVALENT_MACE_DAYS days of index = " r(N)
save "mace_analytic_full.dta", replace


/*==============================================================================
  4. INCIDENCE AND UNADJUSTED KAPLAN-MEIER
==============================================================================*/

cd "$longfmt"
use "mace_analytic_full.dta", clear
drop if prevalent_mace == 1
flowcount "INCIDENT MACE ANALYSIS COHORT"

stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

sts graph, ///
    risktable(0(1)10, size(vsmall) title("Number at risk", size(vsmall))) ///
    xtitle("Years from type 2 diabetes diagnosis", margin(medium)) ///
    ytitle("MACE-5-free survival (%)", margin(r+6)) ///
    xlabel(0(1)10) ///
    ylabel(0.0 "0" 0.2 "20" 0.4 "40" 0.6 "60" 0.8 "80" 1.0 "100") ///
    title("") ///
    graphregion(color(white)) plotregion(color(white)) scheme(s2color)
graph export "mace_km_overall.png", replace width(1200) height(900)

di as txt _n "Overall incidence:"
stptime, per(1000)
di as txt _n "Cumulative incidence:"
sts list, at(5 9.99) failure

*--- component counts and rates on a common denominator ---*
di as txt _n "Component counts among incident MACE:"
foreach comp in acs acs_mi acs_ua stroke hf tia cv_death_first {
    quietly count if failure == 1 & `comp' == 1
    di as txt "  `comp': " r(N)
}
di as txt _n "Component rates per 1,000 person-years:"
foreach comp in acs stroke hf tia cv_death_first {
    quietly stset exit_date, id(patid) origin(time diag_date) ///
        failure(`comp' == 1) scale(365.25)
    quietly stptime, per(1000)
    di as txt "  `comp': " %6.3f r(rate)
}
quietly stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

save "mace_incident_analysis.dta", replace


/*==============================================================================
  5. PRIMARY COMPLETE-CASE COX
==============================================================================*/

cd "$longfmt"
use "mace_incident_analysis.dta", clear
drop if gender == 3
stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

di as txt _n "Cumulative MACE-5 incidence by prescribing status:"
sts list, by(ever_any_drug) at(5 9.99) failure

di as txt _n "Unadjusted:"
stcox i.ever_any_drug, nolog

di as txt _n "Adjusted (primary):"
stcox i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
estimates store cox_primary_cc
estimates save "cox_mace_primary.ster", replace
estat phtest, detail


/*==============================================================================
  6. MACE BY TRAJECTORY
==============================================================================*/

cd "$longfmt"
foreach t in hba1c_trajectory bmi_trajectory {

    local lbl = cond("`t'" == "hba1c_trajectory", "HbA1c", "BMI")
    di as txt _n "{hline 60}"
    di as txt "MACE by `lbl' trajectory"
    di as txt "{hline 60}"

    use "mace_incident_analysis.dta", clear
    drop if gender == 3
    keep if !missing(`t')
    di as txt "  N: " _N

    stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

    tab failure `t', col chi2
    stptime, by(`t') per(1000)
    sts list, at(5 9.99) by(`t') failure
    sts test `t'

    sts graph, by(`t') ///
        xtitle("Years from type 2 diabetes diagnosis", size(medium)) ///
        ytitle("MACE-free survival probability", size(medium)) ///
        legend(order(1 "Improving" 2 "Stable" 3 "Worsening") rows(1) position(6)) ///
        xlabel(0(2)10) ylabel(0.65(0.05)1.0, angle(horizontal)) ///
        risktable title("") ///
        graphregion(color(white)) plotregion(color(white)) scheme(s2color)
    graph export "mace_km_`t'.png", replace width(1200) height(900)

    stcox i.`t' i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
    estimates store cox_`t'
    estimates save "cox_mace_`t'.ster", replace

    stcox i.`t'##i.ever_any_drug i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
    estimates store cox_`t'_int
}

*--- combined model ---*
di as txt _n "Combined HbA1c + BMI trajectory model"
use "mace_incident_analysis.dta", clear
drop if gender == 3
keep if !missing(hba1c_trajectory) & !missing(bmi_trajectory)
stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)
stcox i.hba1c_trajectory i.bmi_trajectory i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
estimates store cox_combined
estat phtest, detail

*--- continuous slopes ---*
di as txt _n "Continuous slope dose-response"
use "mace_incident_analysis.dta", clear
drop if gender == 3
stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

foreach s in hba1c_b_slope bmi_b_slope {
    capture confirm variable `s'
    if _rc continue
    di as txt _n "  `s':"
    summarize `s', detail
    stcox c.`s' i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
    estimates store cox_`s'
}
di as txt _n "  Non-linearity check (quadratic HbA1c slope):"
stcox c.hba1c_b_slope##c.hba1c_b_slope i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron

di as txt _n "  Both slopes:"
stcox c.hba1c_b_slope c.bmi_b_slope i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
estimates store cox_bothslopes


/*==============================================================================
  7. COMPETING-RISK CUMULATIVE INCIDENCE

  Non-cardiovascular deaths count as competing events only if they occur
  within follow-up (i.e. on or before exit_date); deaths recorded after the
  follow-up cap are not treated as competing events, which would otherwise
  deflate the cumulative incidence.
==============================================================================*/

cd "$longfmt"
foreach t in hba1c_trajectory bmi_trajectory {

    use "mace_incident_analysis.dta", clear
    drop if gender == 3
    keep if !missing(`t')

    gen byte cr_event = 0
    replace cr_event = 1 if failure == 1
    replace cr_event = 2 if failure == 0 & !missing(death_date) & death_date <= exit_date

    stset exit_date, id(patid) origin(time diag_date) failure(cr_event == 1) scale(365.25)
    stcrreg i.`t' i.ever_any_drug $COVARS c.c_hba c.c_bmi, compete(cr_event == 2)
    estimates store crr_`t'

    stcurve, cif at(`t' = (1 2 3)) outfile("cif_`t'", replace) ///
        ytitle("Cumulative MACE incidence") ///
        xtitle("Years from type 2 diabetes diagnosis") ///
        legend(order(1 "Improving" 2 "Stable" 3 "Worsening"))

    preserve
        use "cif_`t'.dta", clear
        foreach y in 5 10 {
            gen d`y' = abs(_t - `y')
            sort d`y'
            di as txt "  adjusted CIF at `y' years:"
            list _t ci1 ci2 ci3 in 1, clean noobs
        }
    restore
}


/*==============================================================================
  8. MISSING DATA CHARACTERISATION
==============================================================================*/

di as txt _n "{hline 78}"
di as txt "8: Missing data characterisation"
di as txt "{hline 78}"

cd "$longfmt"
use "mace_incident_analysis.dta", clear

gen byte completeness = .
replace completeness = 1 if !missing(baseline_hba1c) & !missing(baseline_bmi)
replace completeness = 2 if !missing(baseline_hba1c) &  missing(baseline_bmi)
replace completeness = 3 if  missing(baseline_hba1c) & !missing(baseline_bmi)
replace completeness = 4 if  missing(baseline_hba1c) &  missing(baseline_bmi)
label define complab 1 "Both observed" 2 "HbA1c only" 3 "BMI only" 4 "Neither", replace
label values completeness complab
tab completeness

gen byte has_trajectory = !missing(hba1c_trajectory)
tab has_trajectory completeness, col

di as txt _n "Age by completeness group (ANOVA):"
oneway age completeness, tabulate
di as txt _n "Engagement by completeness group (Kruskal-Wallis):"
kwallis engagement, by(completeness)
foreach v in gender smoking_status imd ever_any_drug ihd_prev hf_prev af_prev {
    di as txt _n "`v' by completeness group:"
    tab `v' completeness, col chi2
}

di as txt _n "Multinomial model for completeness group:"
mlogit completeness i.gender c.age i.smoking_status i.imd c.engagement ///
    i.ever_any_drug i.ihd_prev i.hf_prev i.af_prev, baseoutcome(1) rrr nolog


/*==============================================================================
  9. MULTIPLE IMPUTATION
==============================================================================*/

di as txt _n "{hline 78}"
di as txt "9: Multiple imputation ($MI_N datasets)"
di as txt "{hline 78}"

cd "$longfmt"
use "mace_incident_analysis.dta", clear
drop if gender == 3

*--- predictors registered as complete must actually be complete ---*
local badvars ""
foreach v in age gender smoking_status imd ihd_prev hf_prev af_prev ///
             htn_prev copd_prev ckd_prev highchol_prev ever_any_drug ///
             failure diag_date exit_date {
    quietly count if missing(`v')
    if r(N) > 0 local badvars "`badvars' `v'"
}
if "`badvars'" != "" {
    di as error "STOP: registered-complete variables contain missing values:`badvars'"
    di as error "mi impute would drop these observations silently."
    exit 459
}

mi set wide
mi register imputed baseline_hba1c baseline_bmi
mi register regular age gender smoking_status imd ihd_prev hf_prev af_prev ///
    htn_prev copd_prev ckd_prev highchol_prev ever_any_drug ///
    failure diag_date exit_date

mi impute chained ///
    (regress) baseline_hba1c ///
    (regress) baseline_bmi   ///
    = i.gender c.age i.smoking_status i.imd ///
      ihd_prev hf_prev af_prev htn_prev copd_prev ckd_prev highchol_prev ///
      i.ever_any_drug failure c.exit_date ///
    , add($MI_N) rseed($MI_SEED) burnin($MI_BURNIN) dots force

quietly mi xeq 0: summarize baseline_hba1c
local mean_hba1c = r(mean)
quietly mi xeq 0: summarize baseline_bmi
local mean_bmi = r(mean)
mi passive: gen c_hba_mi = baseline_hba1c - `mean_hba1c'
mi passive: gen c_bmi_mi = baseline_bmi   - `mean_bmi'

save "mace_mi_imputed.dta", replace

mi stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)
mi estimate, hr dots: stcox i.ever_any_drug $COVARS c.c_hba_mi c.c_bmi_mi
estimates store cox_primary_mi
estimates save "cox_mace_primary_mi.ster", replace


/*==============================================================================
  10. TRAJECTORY THRESHOLD SENSITIVITY
==============================================================================*/

di as txt _n "{hline 78}"
di as txt "10: Trajectory threshold sensitivity"
di as txt "{hline 78}"

cd "$longfmt"
foreach s in hba1c bmi {

    local slope = cond("`s'" == "hba1c", "hba1c_b_slope", "bmi_b_slope")
    local cuts  = cond("`s'" == "hba1c", "0 0.5 1 2 3 5", "-1 0 0.5 1 2 3")
    local unit  = cond("`s'" == "hba1c", "mmol/mol/yr", "kg/m2/yr")

    di as txt _n "--- `s' thresholds ---"
    matrix `s'_thr = J(6, 6, .)
    local row = 1

    foreach t of local cuts {
        use "mace_incident_analysis.dta", clear
        drop if gender == 3
        keep if !missing(`slope')

        gen byte worsen = (`slope' > `t')
        quietly count if worsen == 1
        local nw = r(N)

        stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)
        quietly stcox i.worsen i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron

        local hr  = exp(_b[1.worsen])
        local lci = exp(_b[1.worsen] - 1.96 * _se[1.worsen])
        local uci = exp(_b[1.worsen] + 1.96 * _se[1.worsen])
        local p   = 2 * (1 - normal(abs(_b[1.worsen] / _se[1.worsen])))

        di as txt "  slope > `t' `unit': N=`nw', HR=" %5.3f `hr' ///
            " (" %5.3f `lci' "-" %5.3f `uci' "), p=" %6.4f `p'

        matrix `s'_thr[`row',1] = `t'
        matrix `s'_thr[`row',2] = `nw'
        matrix `s'_thr[`row',3] = `hr'
        matrix `s'_thr[`row',4] = `lci'
        matrix `s'_thr[`row',5] = `uci'
        matrix `s'_thr[`row',6] = `p'
        local row = `row' + 1
    }

    matrix colnames `s'_thr = Threshold N_worsening HR LCI UCI p
    matrix list `s'_thr, format(%9.3f)

    preserve
        clear
        svmat `s'_thr, names(col)
        export delimited using "threshold_sensitivity_`s'.csv", replace
    restore
}


/*==============================================================================
  11. TIME-VARYING PRESCRIBING SENSITIVITY

  Ever/never classification of a time-dependent exposure carries immortal time
  bias: a patient must survive event-free long enough to receive a first
  prescription in order to be classified as exposed.
==============================================================================*/

di as txt _n "{hline 78}"
di as txt "11: Time-varying prescribing sensitivity"
di as txt "{hline 78}"

cd "$longfmt"
use "mace_incident_analysis.dta", clear
drop if gender == 3
keep if !missing(hba1c_trajectory)

* prescriptions before the index date count as exposed from t=0
gen double rx_start = max(first_drug_date, diag_date) if ever_any_drug == 1
format rx_start %td

stset exit_date, id(patid) origin(time diag_date) failure(failure == 1) scale(365.25)

di as txt _n "(a) as published, fixed exposure:"
stcox i.hba1c_trajectory i.ever_any_drug $COVARS c.c_hba c.c_bmi, nolog efron
estimates store traj_fixed

stsplit period, after(time = rx_start) at(0)
gen byte on_drug = (period == 0)
* stsplit leaves period==0 on unsplit records, so never-treated patients would
* otherwise be coded as exposed for their whole follow-up
replace on_drug = 0 if ever_any_drug == 0

di as txt _n "Verification (both must be clean):"
tab on_drug if ever_any_drug == 0
bysort patid: egen byte ever_on = max(on_drug)
tab ever_on ever_any_drug

di as txt _n "(b) time-varying exposure:"
stcox i.hba1c_trajectory i.on_drug $COVARS c.c_hba c.c_bmi, nolog efron
estimates store traj_tv
estimates save "cox_mace_traj_tv.ster", replace

estimates table traj_fixed traj_tv, b(%7.4f) se eform


log close

di as txt _n "06_mace complete."
