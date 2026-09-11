/*==============================================================================
  03_trajectories.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Longitudinal HbA1c and BMI trajectories.

    Section 1  Build the longitudinal analysis file (shared by all models)
    Section 2  Population-averaged trajectories
    Section 3  Patient-level classification via BLUPs
    Section 4  Annual rate of change by subgroup

  INPUT   cohort_long.dta   
          baseline_final.dta                     
  OUTPUT  hba1c_analysis_ready.dta, bmi_analysis_ready.dta
          hba1c_trajectories.dta, bmi_trajectories.dta
          yearly estimate files and figures
==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

global TRAJ_ADJ_ENGAGEMENT 0

log using "$logdir/03_trajectories.log", replace text

di as txt _n "{hline 78}"
di as txt "STEP 03: HbA1c and BMI trajectories"
di as txt "{hline 78}"

cd "$longfmt"


/*==============================================================================
  SECTION 1 — BUILD THE LONGITUDINAL ANALYSIS FILE

  One file per outcome. Both carry baseline HbA1c and BMI as covariates, so
  patients missing either are excluded from both, keeping the two models on a
  common population.
==============================================================================*/

foreach out in hba1c bmi {

    di as txt _n "--- Building analysis file: `out' ---"

    use "cohort_long.dta", clear

    keep if !missing(`out') & !missing(eventdate) & !missing(baseline)
    if "`out'" == "hba1c" keep if inrange(hba1c, $HBA1C_MIN, $HBA1C_MAX)
    if "`out'" == "bmi"   keep if inrange(bmi,   $BMI_MIN,   $BMI_MAX)

    gen t_years = (eventdate - baseline) / 365.25
    drop if t_years < 0
    drop if t_years > $TRAJ_MAXYRS

    * baseline values and covariates
    merge m:1 patid using "baseline_final.dta", ///
        keepusing(baseline_hba1c baseline_bmi gender age smoking_status imd ///
                  ethnicity engagement hes_linked ///
                  htn_prev ihd_prev copd_prev highchol_prev ///
                  af_prev hf_prev ckd_prev) ///
        keep(master match) nogen

    drop if missing(baseline_hba1c) | missing(baseline_bmi)

    * centred baseline covariates
    quietly summarize baseline_hba1c
    gen c_bbha = baseline_hba1c - r(mean)
    quietly summarize baseline_bmi
    gen c_bbmi = baseline_bmi - r(mean)

    * require >=2 observations so a random slope is estimable
    bysort patid: gen int _nobs = _N
    drop if _nobs < 2
    drop _nobs

    quietly unique patid
    di as txt "  observations: " _N "  patients: " r(unique)

    save "`out'_analysis_ready.dta", replace
}


/*==============================================================================
  SECTION 2 — POPULATION-AVERAGED TRAJECTORIES
==============================================================================*/

foreach out in hba1c bmi {

    local unit = cond("`out'" == "hba1c", "mmol/mol", "kg/m2")
    local ylab = cond("`out'" == "hba1c", "HbA1c (mmol/mol)", "BMI (kg/m{sup:2})")

    di as txt _n "{hline 60}"
    di as txt "Population-averaged `out' trajectory"
    di as txt "{hline 60}"

    use "`out'_analysis_ready.dta", clear

    mixed `out' c.t_years c.c_bbha c.c_bbmi ///
        i.gender c.age i.smoking_status i.imd c.engagement ///
        i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
        i.af_prev i.hf_prev i.ckd_prev ///
        || patid: t_years, covariance(unstructured) reml

    estimates store `out'_popavg
    estimates save "`out'_popavg.ster", replace

    *--- adjusted marginal means at each year ---*
    margins, at(t_years=(0(1)$TRAJ_MAXYRS))
    matrix M = r(table)'

    preserve
        clear
        svmat M, names(col)
        gen year = _n - 1
        keep year b ll ul
        rename (b ll ul) (mean lci uci)
        save "`out'_yearly_estimates.dta", replace
        list, clean noobs
    restore

    *--- trajectory figure ---*
    preserve
        use "`out'_yearly_estimates.dta", clear
        twoway (rarea lci uci year, color(gs12)) ///
               (line mean year, lcolor(navy) lwidth(medthick)), ///
            xtitle("Years from type 2 diabetes diagnosis") ///
            ytitle("`ylab'") ///
            xlabel(0(1)$TRAJ_MAXYRS) ///
            legend(off) ///
            graphregion(color(white)) plotregion(color(white))
        graph export "`out'_trajectory.png", replace width(1200) height(800)
    restore
}


/*==============================================================================
  SECTION 3 — PATIENT-LEVEL CLASSIFICATION

  The random slope is the patient's rate of change RELATIVE TO the cohort
  average, not their absolute rate.
==============================================================================*/

foreach out in hba1c bmi {

    local unit = cond("`out'" == "hba1c", "mmol/mol/yr", "kg/m2/yr")

    di as txt _n "{hline 60}"
    di as txt "Patient-level `out' trajectory classification"
    di as txt "{hline 60}"

    use "`out'_analysis_ready.dta", clear

    mixed `out' c.t_years c.c_bbha c.c_bbmi ///
        i.gender c.age i.smoking_status i.imd c.engagement ///
        i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev ///
        i.af_prev i.hf_prev i.ckd_prev ///
        || patid: t_years, covariance(unstructured) reml

    predict b_intercept b_slope, reffects
    predict b_intercept_se b_slope_se, reses

    * fixed-effect slope, so the absolute rate can also be reported
    scalar fe_slope = _b[t_years]
    di as txt "  cohort-average slope: " %6.3f fe_slope " `unit'"

    bysort patid: keep if _n == 1

    *--- classification: relative to the cohort average -----------------------
    gen byte trajectory = .
    replace trajectory = 1 if b_slope < -$TRAJ_CUT             & !missing(b_slope)
    replace trajectory = 2 if inrange(b_slope, -$TRAJ_CUT, $TRAJ_CUT) & !missing(b_slope)
    replace trajectory = 3 if b_slope > $TRAJ_CUT              & !missing(b_slope)

    *--- units-correct labels, defined separately per outcome -----------------
    if "`out'" == "hba1c" {
        label define traj_hba1c ///
            1 "Improving (>1 mmol/mol/yr below cohort average)" ///
            2 "Stable (within 1 mmol/mol/yr of cohort average)" ///
            3 "Worsening (>1 mmol/mol/yr above cohort average)", replace
        label values trajectory traj_hba1c
    }
    else {
        label define traj_bmi ///
            1 "Improving (>1 kg/m2/yr below cohort average)" ///
            2 "Stable (within 1 kg/m2/yr of cohort average)" ///
            3 "Worsening (>1 kg/m2/yr above cohort average)", replace
        label values trajectory traj_bmi
    }

    *--- absolute individual slope, for the dose-response analyses ------------
    gen indiv_slope = fe_slope + b_slope
    label variable indiv_slope "Absolute individual slope (`unit')"

    quietly count if missing(b_slope)
    di as txt "  unclassified (no estimable slope): " r(N)
    tab trajectory, missing

    keep patid trajectory b_intercept b_slope b_intercept_se b_slope_se indiv_slope
    save "`out'_trajectories.dta", replace
}


/*==============================================================================
  SECTION 4 — ANNUAL RATE OF CHANGE BY SUBGROUP

  Restricted to 0-$TRAJ_MAXYRS years
==============================================================================*/

*--- covariate string, with optional engagement term ---*
local RATE_COVARS "i.gender c.age i.smoking_status i.imd i.htn_prev i.ihd_prev i.copd_prev i.highchol_prev i.ckd_prev i.af_prev i.hf_prev"
if $TRAJ_ADJ_ENGAGEMENT == 1 local RATE_COVARS "`RATE_COVARS' c.engagement"

foreach out in hba1c bmi {

    local unit = cond("`out'" == "hba1c", "mmol/mol/year", "kg/m2/year")
    local ylab = cond("`out'" == "hba1c", "HbA1c (mmol/mol)", "BMI (kg/m{sup:2})")

    di as txt _n "{hline 78}"
    di as txt "Annual rate of change in `out' by subgroup"
    di as txt "{hline 78}"

    use "`out'_analysis_ready.dta", clear

    *--- subgroup definitions ---*
    gen byte age_band = .
    replace age_band = 1 if age < 40
    replace age_band = 2 if inrange(age, 40, 49)
    replace age_band = 3 if inrange(age, 50, 59)
    replace age_band = 4 if inrange(age, 60, 69)
    replace age_band = 5 if age >= 70 & !missing(age)
    label define ageband 1 "<40" 2 "40-49" 3 "50-59" 4 "60-69" 5 ">=70", replace
    label values age_band ageband

    gen byte imd_tertile = .
    replace imd_tertile = 1 if inrange(imd, 1, 3)
    replace imd_tertile = 2 if inrange(imd, 4, 7)
    replace imd_tertile = 3 if inrange(imd, 8, 10)
    label define imdt 1 "Least deprived (1-3)" 2 "Middle (4-7)" 3 "Most deprived (8-10)", replace
    label values imd_tertile imdt

    /*--------------------------------------------------------------------
      Overall slope
    --------------------------------------------------------------------*/
    di as txt _n "--- `out': overall ---"
    mixed `out' c.t_years `RATE_COVARS' ///
        || patid: t_years, covariance(unstructured) reml
    lincom t_years
    di as txt "Overall `out' slope: " %6.3f r(estimate) ///
        " (95% CI " %6.3f r(lb) " to " %6.3f r(ub) ") `unit'"
    estimates store `out'_overall

    /*--------------------------------------------------------------------
      Subgroup slopes
      Ethnicity excludes the "Unknown" category (6) from the model.
    --------------------------------------------------------------------*/
    foreach sg in age_band gender imd_tertile ethnicity {

        di as txt _n "--- `out': by `sg' ---"

        preserve
            if "`sg'" == "ethnicity" drop if ethnicity == 6

            quietly levelsof `sg', local(levels)
            local ref : word 1 of `levels'

            *--- drop the subgroup's own term from the adjustment set
            local adj "`RATE_COVARS'"
            if "`sg'" == "age_band"    local adj : subinstr local adj "c.age" "", all
            if "`sg'" == "gender"      local adj : subinstr local adj "i.gender" "", all
            if "`sg'" == "imd_tertile" local adj : subinstr local adj "i.imd" "", all

            mixed `out' c.t_years c.t_years#i.`sg' i.`sg' `adj' ///
                || patid: t_years, covariance(unstructured) reml
            estimates store `out'_`sg'

            *--- slope within each level ---*
            di as txt "Slopes by `sg' (`unit'):"
            foreach l of local levels {
                if `l' == `ref' quietly lincom t_years
                else            quietly lincom t_years + `l'.`sg'#c.t_years
                di as txt "  level `l': " %7.3f r(estimate) ///
                    " (95% CI " %7.3f r(lb) " to " %7.3f r(ub) ")"
            }

            *--- joint test of slope heterogeneity ---*
            testparm c.t_years#i.`sg'
            di as txt "  joint test of slope differences: p = " %6.4f r(p)

            *--- predicted trajectories ---*
            margins `sg', at(t_years=(0(1)$TRAJ_MAXYRS))
            marginsplot, ///
                xtitle("Years from type 2 diabetes diagnosis") ///
                ytitle("`ylab'") ///
                title("") ///
                recast(line) recastci(rarea) ciopts(color(%20)) ///
                xlabel(0(1)$TRAJ_MAXYRS) ///
                legend(position(6) rows(1)) ///
                graphregion(color(white)) plotregion(color(white))
            graph export "`out'_rate_of_change_by_`sg'.png", replace width(1200) height(800)
        restore
    }

    save "`out'_rate_of_change.dta", replace
}

log close

di as txt _n "03_trajectories complete."
di as txt "  Patient-level groups: hba1c_trajectories.dta / bmi_trajectories.dta"
di as txt "  Check the unclassified counts reported in Section 3 against the"
di as txt "  numbers reported in the manuscript."
