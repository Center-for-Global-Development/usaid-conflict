/*
01_build_panel.do
Builds country-month panel merging ACLED conflict data with GODAD USAID disbursements.
Output: $output/acled_conflict_panel.dta
*/

clear all
set more off

// SET THIS TO YOUR PROJECT ROOT
global root "/Users/lee/Library/CloudStorage/Dropbox-CGDEducation/Lee Crawfurd/blogs/2026-02-usaid-conflict"

global input   "$root/data/input"
global output  "$root/data"
global figures "$root/figures"
cap mkdir "$output"
cap mkdir "$figures"

// ============================================================
// 1. ACLED Africa data
// ============================================================

import delimited "$input/Africa_ACLED.csv", clear varnames(1)

// Parse date
gen date = date(week, "YMD")
format date %td
gen year  = year(date)
gen month = month(date)
gen ym    = ym(year, month)
format ym %tm

// Conflict category flags
gen is_armed   = inlist(event_type, "Battles", "Explosions/Remote violence")
gen is_unrest  = inlist(event_type, "Protests", "Riots")
gen is_vac     = (event_type == "Violence against civilians")
gen is_militia = inlist(sub_event_type, "Mob violence", "Looting/property destruction")

// Aggregate to country-month: all events
preserve
	collapse (sum) events_all=events fat_all=fatalities, by(country year month ym)
	tempfile cm_all
	save `cm_all'
restore

// Armed conflict
preserve
	keep if is_armed == 1
	collapse (sum) events_armed=events fat_armed=fatalities, by(country year month ym)
	tempfile cm_armed
	save `cm_armed'
restore

// Social unrest
preserve
	keep if is_unrest == 1
	collapse (sum) events_unrest=events fat_unrest=fatalities, by(country year month ym)
	tempfile cm_unrest
	save `cm_unrest'
restore

// Violence against civilians
preserve
	keep if is_vac == 1
	collapse (sum) events_vac=events fat_vac=fatalities, by(country year month ym)
	tempfile cm_vac
	save `cm_vac'
restore

// Militia/mob violence
preserve
	keep if is_militia == 1
	collapse (sum) events_militia=events fat_militia=fatalities, by(country year month ym)
	tempfile cm_militia
	save `cm_militia'
restore

// Merge all categories
use `cm_all', clear
merge 1:1 country ym using `cm_armed',  nogen
merge 1:1 country ym using `cm_unrest', nogen
merge 1:1 country ym using `cm_vac',    nogen
merge 1:1 country ym using `cm_militia', nogen

// Fill missing with zero
foreach v of varlist events_* fat_* {
	replace `v' = 0 if `v' == .
}

// ============================================================
// 2. USAID treatment from GODAD
// ============================================================

preserve
	use "$input/Crawfurd_adm1.dta", clear

	// Aggregate to country-year
	collapse (sum) USA_disb, by(name_0 year)

	// Average 2015-2020
	keep if year >= 2015 & year <= 2020
	collapse (mean) usaid_avg = USA_disb, by(name_0)

	// Harmonise country names to match ACLED
	gen country = name_0
	replace country = "Democratic Republic of Congo" if name_0 == "Democratic Republic of the Congo"
	replace country = "Ivory Coast" if name_0 == "Côte d'Ivoire"
	replace country = "eSwatini" if name_0 == "Swaziland"
	replace country = "Sao Tome and Principe" if name_0 == "São Tomé and Príncipe"

	keep country usaid_avg
	tempfile treatment
	save `treatment'
restore

merge m:1 country using `treatment', nogen
replace usaid_avg = 0 if usaid_avg == .

// ============================================================
// 3. Population (approx 2020, millions)
// ============================================================

gen pop = .
replace pop = 206   if country == "Nigeria"
replace pop = 115   if country == "Ethiopia"
replace pop = 90    if country == "Democratic Republic of Congo"
replace pop = 60    if country == "Tanzania"
replace pop = 59    if country == "South Africa"
replace pop = 54    if country == "Kenya"
replace pop = 46    if country == "Uganda"
replace pop = 44    if country == "Sudan"
replace pop = 44    if country == "Algeria"
replace pop = 37    if country == "Morocco"
replace pop = 33    if country == "Angola"
replace pop = 31    if country == "Mozambique"
replace pop = 31    if country == "Ghana"
replace pop = 28    if country == "Madagascar"
replace pop = 27    if country == "Cameroon"
replace pop = 26    if country == "Ivory Coast"
replace pop = 24    if country == "Niger"
replace pop = 21    if country == "Burkina Faso"
replace pop = 20    if country == "Mali"
replace pop = 19    if country == "Malawi"
replace pop = 18    if country == "Zambia"
replace pop = 17    if country == "Senegal"
replace pop = 16    if country == "Chad"
replace pop = 16    if country == "Somalia"
replace pop = 15    if country == "Zimbabwe"
replace pop = 13    if country == "Guinea"
replace pop = 13    if country == "Rwanda"
replace pop = 12    if country == "Benin"
replace pop = 12    if country == "Burundi"
replace pop = 12    if country == "Tunisia"
replace pop = 11    if country == "South Sudan"
replace pop = 8     if country == "Togo"
replace pop = 8     if country == "Sierra Leone"
replace pop = 7     if country == "Libya"
replace pop = 5.5   if country == "Congo"
replace pop = 5     if country == "Liberia"
replace pop = 5     if country == "Central African Republic"
replace pop = 5     if country == "Mauritania"
replace pop = 3.5   if country == "Eritrea"
replace pop = 2.4   if country == "Botswana"
replace pop = 2.4   if country == "Gambia"
replace pop = 2.2   if country == "Gabon"
replace pop = 2.1   if country == "Lesotho"
replace pop = 2     if country == "Guinea-Bissau"
replace pop = 1.4   if country == "Equatorial Guinea"
replace pop = 1.3   if country == "Mauritius"
replace pop = 1.2   if country == "eSwatini"
replace pop = 1     if country == "Djibouti"
replace pop = 0.87  if country == "Comoros"
replace pop = 0.56  if country == "Cape Verde"
replace pop = 0.1   if country == "Seychelles"
replace pop = 102   if country == "Egypt"

drop if pop == .

// ============================================================
// 4. Construct analysis variables
// ============================================================

// Events per million + IHS
foreach cat in all armed unrest vac militia {
	gen epm_`cat' = events_`cat' / pop
	gen ihs_`cat' = asinh(epm_`cat')
}

// USAID per capita (USD per person) and standardised
gen usaid_pc = usaid_avg / (pop * 1000000)
sum usaid_pc
gen usaid_pc_std = usaid_pc / r(sd)

// Event time: Jan 2025 = month 0
gen event_time = (year - 2025) * 12 + month - 1

// Post indicator
gen post = (event_time >= 0)

// Country numeric ID for FE
encode country, gen(cid)

// Label key variables
label var events_all     "Total conflict events"
label var events_armed   "Armed conflict events (battles + explosions)"
label var events_unrest  "Social unrest events (protests + riots)"
label var events_vac     "Violence against civilians"
label var events_militia "Militia/mob violence events"
label var usaid_pc       "USAID per capita (USD)"
label var usaid_pc_std   "USAID per capita (standardised)"
label var event_time     "Months relative to USAID suspension (Jan 2025 = 0)"
label var post           "Post-suspension indicator"

save "$output/acled_conflict_panel.dta", replace

di "Panel built: " _N " observations, " r(N) " countries"
tab country if event_time == 0, m
