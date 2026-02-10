/*
05_country_breakdown.do
Scatterplots of Δ conflict vs USAID exposure.
- Subnational (admin1): Δ armed fatalities vs USAID, Δ armed events vs USAID
- Country-level: same
Requires: subnational_panel.dta, acled_conflict_panel.dta
*/

clear all
set more off

// SET THIS TO YOUR PROJECT ROOT
global root "/Users/lee/Library/CloudStorage/Dropbox-CGDEducation/Lee Crawfurd/blogs/2026-02-usaid-conflict"

global input   "$root/data/input"
global output  "$root/data"
global figures "$root/figures"

// ============================================================
// PART A: SUBNATIONAL SCATTERPLOTS
// ============================================================
use "$input/subnational_panel.dta", clear

// Pre-period: event_time -6 to -1
preserve
keep if event_time >= -6 & event_time <= -1
collapse (mean) events_armed_pre = events_armed fat_armed_pre = fat_armed, ///
	by(region_id COUNTRY ADMIN1 usaid_avg)
tempfile pre
save `pre'
restore

// Post-period: event_time 0 to 12
preserve
keep if event_time >= 0 & event_time <= 12
collapse (mean) events_armed_post = events_armed fat_armed_post = fat_armed, ///
	by(region_id)
tempfile post_data
save `post_data'
restore

use `pre', clear
merge 1:1 region_id using `post_data', nogen

gen event_change = events_armed_post - events_armed_pre
gen fat_change = fat_armed_post - fat_armed_pre

// Log USAID for x-axis (many zeros, so use asinh)
gen ihs_usaid = asinh(usaid_avg)

// Label for large outliers
gen label = ""
// Find top regions by fatality change
gsort -fat_change
replace label = ADMIN1 + " (" + COUNTRY + ")" in 1/5
gsort fat_change
replace label = ADMIN1 + " (" + COUNTRY + ")" in 1/3

// --- Scatter: Δ fatalities vs USAID (subnational) ---
twoway (scatter fat_change usaid_avg if usaid_avg > 0 & label == "", ///
		mcolor(gs12) msymbol(circle) msize(vsmall) jitter(3)) ///
	(scatter fat_change usaid_avg if usaid_avg > 0 & label != "", ///
		mcolor(cranberry) msymbol(circle) msize(small) ///
		mlabel(label) mlabposition(3) mlabsize(tiny) mlabcolor(black)) ///
	(lfit fat_change usaid_avg if usaid_avg > 0, lcolor(cranberry) lwidth(medthick)), ///
	yline(0, lcolor(black) lwidth(thin)) ///
	title("Change in Armed Conflict Fatalities vs USAID Disbursements") ///
	subtitle("Each point is an admin1 region") ///
	ytitle("Δ fatalities/month (2025 vs Jul–Dec 2024)") ///
	xtitle("USAID disbursements (avg 2015–2020, millions USD)") ///
	graphregion(color(white)) plotregion(margin(small)) ///
	legend(off) ///
	note("Admin1 regions with positive USAID shown. Line shows linear fit.")
graph export "$figures/scatter_subnational_fatalities.png", replace width(1200)

// --- Scatter: Δ events vs USAID (subnational) ---
twoway (scatter event_change usaid_avg if usaid_avg > 0 & label == "", ///
		mcolor(gs12) msymbol(circle) msize(vsmall) jitter(3)) ///
	(scatter event_change usaid_avg if usaid_avg > 0 & label != "", ///
		mcolor(cranberry) msymbol(circle) msize(small) ///
		mlabel(label) mlabposition(3) mlabsize(tiny) mlabcolor(black)) ///
	(lfit event_change usaid_avg if usaid_avg > 0, lcolor(cranberry) lwidth(medthick)), ///
	yline(0, lcolor(black) lwidth(thin)) ///
	title("Change in Armed Conflict Events vs USAID Disbursements") ///
	subtitle("Each point is an admin1 region") ///
	ytitle("Δ events/month (2025 vs Jul–Dec 2024)") ///
	xtitle("USAID disbursements (avg 2015–2020, millions USD)") ///
	graphregion(color(white)) plotregion(margin(small)) ///
	legend(off) ///
	note("Admin1 regions with positive USAID shown. Line shows linear fit.")
graph export "$figures/scatter_subnational_events.png", replace width(1200)

// Summary stats
di _newline "============================================"
di "SUBNATIONAL SCATTER SUMMARY"
di "============================================"
su usaid_avg event_change fat_change if usaid_avg > 0
count if usaid_avg > 0
di "Regions with USAID > 0: " r(N)

// Correlation
corr fat_change usaid_avg if usaid_avg > 0
corr event_change usaid_avg if usaid_avg > 0

// ============================================================
// PART B: COUNTRY-LEVEL SCATTERPLOTS (keep for reference)
// ============================================================
use "$output/acled_conflict_panel.dta", clear

// Pre: Jul-Dec 2024
preserve
keep if event_time >= -6 & event_time <= -1
collapse (mean) events_pre = events_armed fat_pre = fat_armed, by(country usaid_pc pop)
tempfile cpre
save `cpre'
restore

// Post: Jan 2025+
preserve
keep if event_time >= 0 & event_time <= 12
collapse (mean) events_post = events_armed fat_post = fat_armed, by(country)
tempfile cpost
save `cpost'
restore

use `cpre', clear
merge 1:1 country using `cpost', nogen
gen event_change = events_post - events_pre
gen fat_change = fat_post - fat_pre

// Export for Python bar chart
export delimited "$output/country_breakdown.csv", replace

// Labels
gen label = ""
replace label = "DRC" if country == "Democratic Republic of Congo"
replace label = "Somalia" if country == "Somalia"
replace label = "S. Sudan" if country == "South Sudan"
replace label = "Mali" if country == "Mali"
replace label = "Nigeria" if country == "Nigeria"
replace label = "Sudan" if country == "Sudan"
replace label = "Ethiopia" if country == "Ethiopia"
replace label = "Kenya" if country == "Kenya"
replace label = "Cameroon" if country == "Cameroon"
replace label = "CAR" if country == "Central African Republic"
replace label = "B. Faso" if country == "Burkina Faso"

twoway (scatter fat_change usaid_pc if label == "", ///
		mcolor(gs10) msymbol(circle) msize(small)) ///
	(scatter fat_change usaid_pc if label != "", ///
		mcolor(cranberry) msymbol(circle) msize(medsmall) ///
		mlabel(label) mlabposition(3) mlabsize(vsmall) mlabcolor(black)) ///
	(lfit fat_change usaid_pc, lcolor(cranberry%50) lpattern(dash)), ///
	yline(0, lcolor(black) lwidth(thin)) ///
	title("Change in Armed Conflict Fatalities vs USAID Exposure") ///
	subtitle("Each point is a country") ///
	ytitle("Δ fatalities/month (2025 vs Jul–Dec 2024)") ///
	xtitle("USAID per capita ($, 2015–2020 avg)") ///
	graphregion(color(white)) plotregion(margin(small)) ///
	legend(off)
graph export "$figures/scatter_fatalities_usaid.png", replace width(1200)

di _newline "Done — all scatterplots."
