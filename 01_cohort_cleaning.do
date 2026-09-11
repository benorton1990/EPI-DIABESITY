/*==============================================================================
  01_cohort_cleaning.do
  <STUDY NAME>  |  CPRD Aurum  |  ISAC/REC protocol <PROTOCOL ID>

  Builds the cleaned long-format analysis file from the raw Aurum extract.

  INPUT   Observation extracts (per subpart), patient file, drug issue file,
          consultation file, medcode lookup
  OUTPUT  part<p>_<sub>_long.dta   (per subpart, in $longfmt)
          cohort_long.dta          (all subparts appended)
==============================================================================*/

if "$root" == "" do "00_config.do"   // run from the repo root, or edit to an absolute path

log using "$logdir/01_cohort_cleaning.log", replace text

di as txt _n "{hline 78}"
di as txt "STEP 01: Cohort cleaning"
di as txt "{hline 78}"

* Which extract parts to process
global PARTS "1"                 // e.g. "1 2 3" if further extracts are added


foreach part of global PARTS {

    local obsdir "$root/raw_extract_part_`part'/observation"
    local files : dir "`obsdir'" files "raw_observation_extract_part`part'_???.dta"

    foreach f of local files {

        local sub = substr("`f'", -7, 3)
        di as txt _n "{hline 78}"
        di as txt "Part `part', subpart `sub'"
        di as txt "{hline 78}"

        /*======================================================================
          1. LOAD OBSERVATIONS AND ATTACH MEDCODE CLASSIFICATION
             medcode_lookup file derived from the Aurum codebook
        ======================================================================*/
        use "`obsdir'/`f'", clear

        merge m:1 medcodeid using "$root/raw_medcode_lookup.dta"
        keep if _merge == 3
        drop _merge
        drop consid parentobsid probobsid

        /*======================================================================
          2. LINK TO PATIENT FILE
        ======================================================================*/
        merge m:1 patid using "$root/raw_patient_extract_part`part'.dta"
        keep if _merge == 3
        drop _merge mob

        if _N == 0 {
            di as txt "No matched patients in `f' - skipping"
            continue
        }
        flowcount "loaded and linked to patient file"

        /*======================================================================
          3. COHORT DEFINITION
        ======================================================================*/

        *--- phenotype flags ---*
        gen byte has_t2d       = category == "t2d"
        gen byte has_obesity   = category == "obesity"
        gen byte has_t1d       = category == "t1d"
        gen byte has_bariatric = category == "bs"

        *--- require T2D and obesity; exclude any T1D code (patient level) ---*
        egen byte patient_has_t2d     = max(has_t2d),     by(patid)
        egen byte patient_has_obesity = max(has_obesity), by(patid)
        egen byte patient_has_t1d     = max(has_t1d),     by(patid)

        drop if patient_has_t2d == 0
        flowcount "with a T2D code"
        drop if patient_has_obesity == 0
        flowcount "  and an obesity code"
        drop if patient_has_t1d == 1
        flowcount "  excluding any T1D code"

        *--- observation date ---*
        gen obsdate_num = date(obsdate, "DMY")
        format obsdate_num %td
        drop if missing(obsdate_num)

        *--- earliest date per phenotype ---*
        egen earliest_t2d       = min(cond(has_t2d == 1,       obsdate_num, .)), by(patid)
        egen earliest_obesity   = min(cond(has_obesity == 1,   obsdate_num, .)), by(patid)
        egen earliest_bariatric = min(cond(has_bariatric == 1, obsdate_num, .)), by(patid)
        format earliest_t2d earliest_obesity earliest_bariatric %td

        *--- ensure a single index date per patient ---*
        egen true_earliest_t2d = min(earliest_t2d), by(patid)
        replace earliest_t2d = true_earliest_t2d
        drop true_earliest_t2d

        *--- bariatric surgery more than $BARIATRIC_WINDOW days before index ---*
        drop if !missing(earliest_bariatric) & ///
                earliest_bariatric < earliest_t2d - $BARIATRIC_WINDOW
        flowcount "  excluding prior bariatric surgery"

        *--- pregnancy within +/- $PREGNANCY_WINDOW days of index ---*
        gen byte has_preg = category == "pregnancy"
        egen earliest_preg_near = min(cond(has_preg == 1 & ///
            inrange(obsdate_num, earliest_t2d - $PREGNANCY_WINDOW, ///
                                 earliest_t2d + $PREGNANCY_WINDOW), ///
            obsdate_num, .)), by(patid)
        drop if !missing(earliest_preg_near)
        flowcount "  excluding peri-diagnosis pregnancy"

        *--- obesity code must not postdate index by > $OBESITY_WINDOW days ---*
        drop if earliest_t2d - earliest_obesity < -$OBESITY_WINDOW
        flowcount "  with obesity recorded before or within 1y of index"

        *--- age at diagnosis ---*
        gen age_at_t2d = year(earliest_t2d) - yob
        drop if age_at_t2d < $MIN_AGE
        flowcount "  aged >=$MIN_AGE at diagnosis"

        *--- registration: >= $MIN_REGISTRATION days before index ---*
        gen regstartdate_num = date(regstartdate, "DMY")
        gen days_from_reg = earliest_t2d - regstartdate_num
        drop if days_from_reg < $MIN_REGISTRATION
        flowcount "  registered >=$MIN_REGISTRATION days pre-index"

        drop regstartdate
        rename regstartdate_num registered
        format registered %td

        gen regenddate_num = date(regenddate, "DMY")
        drop regenddate
        rename regenddate_num deregistered
        format deregistered %td

        *--- date of death ---*
        drop emis_ddate
        rename cprd_ddate death_date
        gen death_date_num = date(death_date, "DMY")
        drop death_date
        rename death_date_num death_date
        format death_date %td

        *--- recruitment window ---*
        drop if earliest_t2d < date("$STUDY_START", "DMY") | ///
                earliest_t2d > date("$STUDY_END", "DMY")
        flowcount "  diagnosed $STUDY_START to $STUDY_END"

        *--- peri-diagnosis solid cancer (+/- $CANCER_WINDOW days) ---*
        gen byte _cancer = 0
        foreach c in breast_ca crc hcc lung_ca og_ca ovarian_ca pan_ca ///
                     prostate_ca rcc uterine_ca {
            replace _cancer = 1 if category == "`c'"
        }
        gen _cancer_date = obsdate_num if _cancer == 1
        bysort patid: egen _earliest_cancer = min(_cancer_date)
        gen byte _exclude_cancer = !missing(_earliest_cancer) ///
            & inrange(_earliest_cancer, earliest_t2d - $CANCER_WINDOW, ///
                                        earliest_t2d + $CANCER_WINDOW)
        drop if _exclude_cancer == 1
        flowcount "  excluding peri-diagnosis solid cancer  << FINAL COHORT"

        drop _cancer _cancer_date _earliest_cancer _exclude_cancer

        save "`obsdir'/part`part'_`sub'_cohort_eligible.dta", replace

        /*======================================================================
          4. BMI AND HbA1c  (HbA1c standardised to mmol/mol)
        ======================================================================*/
        gen bmi = value   if category == "bmi"
        gen hba1c = value if category == "hba1c"

        * DCCT (%) values are 4-15; IFCC (mmol/mol) values are >15
        gen str8 hba1c_unit = ""
        replace hba1c_unit = "percent"  if inrange(hba1c, 4, 15)
        replace hba1c_unit = "mmol/mol" if hba1c > 15 & !missing(hba1c)

        gen hba1c_mmol = .
        replace hba1c_mmol = hba1c if hba1c_unit == "mmol/mol"
        replace hba1c_mmol = (hba1c - 2.15) * 10.929 if hba1c_unit == "percent"
        replace hba1c_mmol = round(hba1c_mmol)
        drop hba1c hba1c_unit
        rename hba1c_mmol hba1c

        * Implausible values removed here so every downstream script inherits
        * the same range (see $HBA1C_MIN/$HBA1C_MAX, $BMI_MIN/$BMI_MAX)
        replace hba1c = . if !inrange(hba1c, $HBA1C_MIN, $HBA1C_MAX)
        replace bmi   = . if !inrange(bmi,   $BMI_MIN,   $BMI_MAX)

        /*======================================================================
          5. COMORBIDITIES  (flag on the diagnosis row)
        ======================================================================*/
        foreach cond in ihd htn dyslipid pvd copd ckd af hf anx_dep nafld {
            gen byte `cond' = .
            replace `cond' = 1 if category == "`cond'"
        }

        /*======================================================================
          6. SMOKING STATUS
             0 = never/non, 1 = current, 2 = ex
             Current: code within 4 years before index. Ex: any code before
             index. Where both, the most recent wins.
        ======================================================================*/
        gen time_from_t2d = obsdate_num - earliest_t2d
        gen byte flag_smoker   = (category == "smoker"   & inrange(time_from_t2d, -1460, 0))
        gen byte flag_exsmoker = (category == "exsmoker" & time_from_t2d <= 0)
        gen smoker_date   = obsdate_num if category == "smoker"   & time_from_t2d <= 0
        gen exsmoker_date = obsdate_num if category == "exsmoker" & time_from_t2d <= 0

        bysort patid: egen latest_smoker_date   = max(smoker_date)
        bysort patid: egen latest_exsmoker_date = max(exsmoker_date)
        bysort patid: egen byte ever_smoker     = max(flag_smoker)
        bysort patid: egen byte ever_exsmoker   = max(flag_exsmoker)

        gen byte smoking_status = 0
        replace smoking_status = 2 if ever_exsmoker == 1
        replace smoking_status = 1 if ever_smoker == 1 & ///
            (missing(latest_exsmoker_date) | latest_smoker_date >= latest_exsmoker_date)
        label define smoklab 0 "Never/non-smoker" 1 "Current smoker" 2 "Ex-smoker", replace
        label values smoking_status smoklab

        drop time_from_t2d flag_smoker flag_exsmoker smoker_date exsmoker_date ///
             latest_smoker_date latest_exsmoker_date ever_smoker ever_exsmoker

        /*======================================================================
          7. DEATH
        ======================================================================*/
        gen byte died = !missing(death_date)


        /*======================================================================
          8. HEALTHCARE ENGAGEMENT
             Unique consultation days in the 365 days before index.
        ======================================================================*/
        preserve
            bysort patid: keep if _n == 1
            keep patid earliest_t2d
            format patid %14.0g
            save "`obsdir'/part`part'_`sub'_baseline_dates.dta", replace
        restore

        preserve
            use "$consdir/raw_consultation_dates.dta", clear
            merge m:1 patid using ///
                "`obsdir'/part`part'_`sub'_baseline_dates.dta", keep(match) nogen
            gen byte pre_baseline = inrange(consdate_num, earliest_t2d - 365, earliest_t2d - 1)
            bysort patid consdate_num: gen byte consult_unique = (_n == 1) if pre_baseline == 1
            bysort patid: egen engagement = total(consult_unique)
            collapse (max) engagement, by(patid)
            save "`obsdir'/part`part'_`sub'_engagement.dta", replace
        restore

        merge m:1 patid using "`obsdir'/part`part'_`sub'_engagement.dta", ///
            keep(master match) nogen
        replace engagement = 0 if missing(engagement)

        /*======================================================================
          9. PRESCRIBING
              Drug issues restricted to this subpart's cohort, de-duplicated,
              and filtered to chronic use (>=2 issues) unless pre-baseline.
        ======================================================================*/
        preserve
            use "$root/raw_drug_issue_extract.dta", clear
            gen issuedate_num = date(issuedate, "DMY")
            format issuedate_num %td
            drop issuedate
            drop if issuedate_num < date("$STUDY_START", "DMY")
            drop if issuedate_num > date("31/12/2024", "DMY")
            keep if inlist(class, "metformin","su","insulin","glp","sglt2", ///
                                  "dppiv","tzd","megalitinides","statin")

            merge m:1 patid using ///
                "`obsdir'/part`part'_`sub'_baseline_dates.dta", keep(match) nogen

            count
            if r(N) == 0 {
                * typed empty file so the append below still works
                clear
                gen long   patid = .
                gen str12  class = ""
                gen double issuedate_num = .
                format issuedate_num %td
                gen byte   pre_baseline_drug = .
                gen byte   baseline_window_drug = .
                gen byte   post_baseline_drug = .
                gen double first_drug_date = .
                format first_drug_date %td
                gen byte   has_dosage = .
                gen str100 dosage_text = ""
                save "`obsdir'/part`part'_`sub'_drugissue_clean.dta", replace
            }
            else {
                duplicates drop patid class issuedate_num, force

                * timing relative to index
                gen byte pre_baseline_drug    = (issuedate_num <  earliest_t2d - 90)
                gen byte baseline_window_drug = inrange(issuedate_num, ///
                                                        earliest_t2d - 90, earliest_t2d - 1)
                gen byte post_baseline_drug   = (issuedate_num >= earliest_t2d)

                * first issue date per class
                bysort patid class (issuedate_num): gen byte index_rx = (_n == 1)
                gen index_rx_date = issuedate_num if index_rx == 1
                bysort patid class: egen first_drug_date = min(index_rx_date)
                format first_drug_date %td
                drop index_rx index_rx_date

                * chronic use filter: keep pre-baseline issues, or >=2 issues
                bysort patid class: gen int total_rx = _N
                gen byte chronic_med = (total_rx >= 2)
                drop if chronic_med == 0 & pre_baseline_drug == 0
                drop total_rx chronic_med

                merge m:1 dosageid using "$root/raw_dosage_lookup.dta", keep(master match) nogen
                gen byte has_dosage = !missing(dosage_text)
                drop dosageid drugrecid
                save "`obsdir'/part`part'_`sub'_drugissue_clean.dta", replace
            }
        restore

        append using "`obsdir'/part`part'_`sub'_drugissue_clean.dta"
        sort patid issuedate_num

        foreach d in statin metformin dppiv glp sglt2 su insulin tzd megalitinides {
            gen byte `d' = (class == "`d'")
        }

        /*======================================================================
          10. FINAL LONG-FORMAT CLEANING
        ======================================================================*/

        *--- bariatric and cancer covariate flags ---*
        gen byte bs = (category == "bs")
        foreach c in breast_ca crc hcc lung_ca og_ca ovarian_ca pan_ca ///
                     prostate_ca rcc uterine_ca {
            gen byte `c' = (category == "`c'")
        }
        gen byte cancer_solid = (breast_ca | crc | hcc | lung_ca | og_ca | ///
                                 ovarian_ca | pan_ca | prostate_ca | rcc | uterine_ca)

        *--- drop raw and intermediate columns ---*
        capture drop obsid obsdate enterdate medcodeid numunitid obstypeid ///
            numrangelow numrangehigh term SnomedCTDescriptionId dup_medcode
        capture drop has_t2d has_obesity patient_has_obesity patient_has_t2d ///
            has_t1d patient_has_t1d has_bariatric has_preg earliest_preg_near ///
            earliest_bariatric earliest_obesity days_from_reg
        capture drop registered deregistered
        capture drop issueid prodcodeid quantunitid drug combined
        capture drop dose_frequency dose_interval choice_of_dose dose_max_average ///
            change_dose dose_duration has_dosage
        capture drop pracid

        drop category

        rename died      death
        rename anx_dep   depression
        rename dyslipid  highchol

        *--- single event date from clinical OR prescribing ---*
        assert missing(obsdate_num) | missing(issuedate_num)
        gen double eventdate = cond(!missing(obsdate_num), obsdate_num, issuedate_num)
        format eventdate %td
        drop obsdate_num issuedate_num

        *--- date variable per flagged condition / drug ---*
        foreach var in bmi hba1c ihd htn highchol bs pvd copd ckd af hf depression ///
                       nafld cancer_solid breast_ca crc hcc lung_ca og_ca ovarian_ca ///
                       pan_ca prostate_ca rcc uterine_ca statin metformin dppiv glp ///
                       sglt2 su insulin tzd megalitinides {
            gen `var'_date = eventdate if `var' == 1
        }
        format *_date %td

        *--- propagate patient-level variables across all rows ---*
        foreach var in gender age_at_t2d yob death smoking_status engagement ///
                       prevalent_mace mace_event {
            bysort patid: egen temp = max(`var')
            replace `var' = temp if missing(`var')
            drop temp
        }
        foreach var in death_date mace5_date earliest_t2d {
            capture {
                bysort patid: egen double temp = min(`var')
                replace `var' = temp if missing(`var')
                format `var' %td
                drop temp
            }
        }

        capture drop mace_hosp_date cv_death_date
        capture drop dosage_text daily_dose dose_number dose_unit quantity duration
        capture drop value

        save "$longfmt/part`part'_`sub'_long.dta", replace
        di as txt "Saved subpart `sub'."
    }
}


/*==============================================================================
  11. APPEND ALL SUBPARTS
==============================================================================*/
di as txt _n "{hline 78}"
di as txt "Appending subparts"
di as txt "{hline 78}"

cd "$longfmt"
local allfiles : dir "$longfmt" files "part*_???_long.dta"

local first = 1
foreach f of local allfiles {
    if `first' == 1 {
        use "`f'", clear
        local first = 0
    }
    else {
        append using "`f'"
    }
}

sort patid eventdate
compress
save "cohort_long.dta", replace

flowcount "FINAL cohort, all subparts appended"
di as txt "Observations: " _N

log close
di as txt _n "01_cohort_cleaning complete."
