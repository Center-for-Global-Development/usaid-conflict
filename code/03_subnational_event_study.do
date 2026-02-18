/*
03_subnational_event_study.do
Subnational (admin1) event study coefplots (re-centred: pre-period mean = 0).
Main spec: TWFE (admin1 + month FE). Robustness: admin1 + country×month FE.
Produces: subnational_armed.png, subnational_fatalities_armed.png
Requires: subnational_panel.dta
*/

clear all
set more off

// SET THIS TO YOUR PROJECT ROOT
global root "/Users/lee/Library/CloudStorage/Dropbox-CGDEducation/Lee Crawfurd/blogs/2026-02-usaid-conflict"

global input   "$root/data/input"
global output  "$root/data"
global figures "$root/figures"

use "$input/subnational_panel.dta", clear

// Create IHS fatality variables
gen ihs_fat_armed = asinh(fat_armed)
gen ihs_fat_all = asinh(fat_all)

keep if event_time >= -6 & event_time <= 12

encode region_id, gen(rid)
egen cym = group(country ym)

// Event-time interactions (omit t = -1)
forvalues t = -6/12 {
	if `t' == -1 continue
	local tlab = cond(`t' < 0, "m" + string(abs(`t')), "p" + string(`t'))
	gen et_`tlab' = (event_time == `t') * usaid_std
}

local etvars et_m6 et_m5 et_m4 et_m3 et_m2 et_p0 et_p1 et_p2 et_p3 et_p4 et_p5 et_p6 et_p7 et_p8 et_p9 et_p10 et_p11 et_p12

// ============================================================
// Program: re-centred coefplot via twoway
// ============================================================
capture program drop recentred_plot
program define recentred_plot
	syntax, depvar(string) etvars(string) color(string) ///
		title(string) subtitle(string) filename(string) fespec(string)

	reghdfe `depvar' `etvars', a(`fespec') cl(rid)

	local N = e(N)
	local Ncl = e(N_clust)

	// Pre-treatment mean (5 estimated + reference = 0, so 6 periods)
	local premean = (_b[et_m6] + _b[et_m5] + _b[et_m4] + _b[et_m3] + _b[et_m2] + 0) / 6

	preserve
	clear
	set obs 19
	gen event_time = _n - 7
	gen coef = .
	gen ci_lo = .
	gen ci_hi = .

	// Reference (t = -1, obs 6)
	replace coef  = 0 - `premean' in 6
	replace ci_lo = 0 - `premean' in 6
	replace ci_hi = 0 - `premean' in 6

	// Pre-treatment (obs 1–5)
	local precoefs et_m6 et_m5 et_m4 et_m3 et_m2
	local row = 1
	foreach v of local precoefs {
		replace coef  = _b[`v'] - `premean' in `row'
		replace ci_lo = _b[`v'] - `premean' - 1.96 * _se[`v'] in `row'
		replace ci_hi = _b[`v'] - `premean' + 1.96 * _se[`v'] in `row'
		local ++row
	}

	// Post-treatment (obs 7–19)
	local postcoefs et_p0 et_p1 et_p2 et_p3 et_p4 et_p5 et_p6 et_p7 et_p8 et_p9 et_p10 et_p11 et_p12
	local row = 7
	foreach v of local postcoefs {
		replace coef  = _b[`v'] - `premean' in `row'
		replace ci_lo = _b[`v'] - `premean' - 1.96 * _se[`v'] in `row'
		replace ci_hi = _b[`v'] - `premean' + 1.96 * _se[`v'] in `row'
		local ++row
	}

	gen sig = (ci_lo > 0 | ci_hi < 0)

	twoway (rcap ci_lo ci_hi event_time, lcolor(`color'%30) lwidth(medthin)) ///
		(scatter coef event_time if sig==1, mcolor(`color') msymbol(circle) msize(small)) ///
		(scatter coef event_time if sig==0, mcolor(`color') msymbol(circle_hollow) msize(small)), ///
		yline(0, lcolor(black) lwidth(thin)) ///
		xline(-0.5, lcolor(gray) lpattern(dash)) ///
		title("`title'") ///
		subtitle("`subtitle'") ///
		ytitle("Coefficient (IHS, re-centred)") ///
		xtitle("Months relative to USAID suspension") ///
		xlabel(-6 "Jul 24" -3 "Oct 24" 0 "Jan 25" 3 "Apr 25" ///
			6 "Jul 25" 9 "Oct 25" 12 "Jan 26") ///
		graphregion(color(white)) plotregion(margin(small)) ///
		legend(off) ///
		note("N = `N', Regions = `Ncl'. Filled = p<0.05. Pre-period re-centred to zero.")
	graph export "$figures/`filename'", replace width(1200)

	restore
end

// ============================================================
// Main spec: TWFE (admin1 + month FE)
// ============================================================

// Armed conflict events
recentred_plot, depvar(ihs_events_armed) etvars(`etvars') color(cranberry) ///
	title("Armed Conflict Events — Subnational Event Study") ///
	subtitle("Admin1 FE + Month FE, SEs clustered at admin1") ///
	filename(subnational_armed.png) fespec(rid ym)

// Armed conflict fatalities
recentred_plot, depvar(ihs_fat_armed) etvars(`etvars') color(cranberry) ///
	title("Armed Conflict Fatalities — Subnational Event Study") ///
	subtitle("Admin1 FE + Month FE, SEs clustered at admin1") ///
	filename(subnational_fatalities_armed.png) fespec(rid ym)

// ============================================================
// Two-model comparison: TWFE vs Country×Month FE — Armed Events
// (matches style of twfe_vs_cym_fatalities.png)
// ============================================================

// --- Model 1: TWFE (admin1 + month FE) ---
reghdfe ihs_events_armed `etvars', a(rid ym) cl(rid)

local premean1 = (_b[et_m6] + _b[et_m5] + _b[et_m4] + _b[et_m3] + _b[et_m2] + 0) / 6

preserve
clear
set obs 19
gen event_time = _n - 7
gen coef1 = .
gen ci_lo1 = .
gen ci_hi1 = .

// Reference (t = -1)
replace coef1  = 0 - `premean1' in 6
replace ci_lo1 = 0 - `premean1' in 6
replace ci_hi1 = 0 - `premean1' in 6

local precoefs et_m6 et_m5 et_m4 et_m3 et_m2
local row = 1
foreach v of local precoefs {
	replace coef1  = _b[`v'] - `premean1' in `row'
	replace ci_lo1 = _b[`v'] - `premean1' - 1.96 * _se[`v'] in `row'
	replace ci_hi1 = _b[`v'] - `premean1' + 1.96 * _se[`v'] in `row'
	local ++row
}

local postcoefs et_p0 et_p1 et_p2 et_p3 et_p4 et_p5 et_p6 et_p7 et_p8 et_p9 et_p10 et_p11 et_p12
local row = 7
foreach v of local postcoefs {
	replace coef1  = _b[`v'] - `premean1' in `row'
	replace ci_lo1 = _b[`v'] - `premean1' - 1.96 * _se[`v'] in `row'
	replace ci_hi1 = _b[`v'] - `premean1' + 1.96 * _se[`v'] in `row'
	local ++row
}

tempfile twfe_results
save `twfe_results'
restore

// --- Model 2: Country×Month FE ---
reghdfe ihs_events_armed `etvars', a(rid cym) cl(rid)

local premean2 = (_b[et_m6] + _b[et_m5] + _b[et_m4] + _b[et_m3] + _b[et_m2] + 0) / 6

preserve
use `twfe_results', clear

gen coef2 = .
gen ci_lo2 = .
gen ci_hi2 = .

// Reference (t = -1)
replace coef2  = 0 - `premean2' in 6
replace ci_lo2 = 0 - `premean2' in 6
replace ci_hi2 = 0 - `premean2' in 6

local precoefs et_m6 et_m5 et_m4 et_m3 et_m2
local row = 1
foreach v of local precoefs {
	replace coef2  = _b[`v'] - `premean2' in `row'
	replace ci_lo2 = _b[`v'] - `premean2' - 1.96 * _se[`v'] in `row'
	replace ci_hi2 = _b[`v'] - `premean2' + 1.96 * _se[`v'] in `row'
	local ++row
}

local postcoefs et_p0 et_p1 et_p2 et_p3 et_p4 et_p5 et_p6 et_p7 et_p8 et_p9 et_p10 et_p11 et_p12
local row = 7
foreach v of local postcoefs {
	replace coef2  = _b[`v'] - `premean2' in `row'
	replace ci_lo2 = _b[`v'] - `premean2' - 1.96 * _se[`v'] in `row'
	replace ci_hi2 = _b[`v'] - `premean2' + 1.96 * _se[`v'] in `row'
	local ++row
}

// Slight x-offset for readability
gen x1 = event_time - 0.15
gen x2 = event_time + 0.15

twoway (rcap ci_lo1 ci_hi1 x1, lcolor(navy%30) lwidth(medthin)) ///
	(scatter coef1 x1, mcolor(navy) msymbol(circle) msize(small)) ///
	(rcap ci_lo2 ci_hi2 x2, lcolor(cranberry%30) lwidth(medthin)) ///
	(scatter coef2 x2, mcolor(cranberry) msymbol(triangle) msize(small)), ///
	yline(0, lcolor(black) lwidth(thin) lpattern(dash)) ///
	xline(-0.5, lcolor(gray) lpattern(dash)) ///
	title("Armed Conflict Events: TWFE vs Country×Month FE") ///
	subtitle("Both use admin1 FE, SEs clustered at admin1") ///
	ytitle("Coefficient (IHS events, re-centred)") ///
	xtitle("Months relative to USAID suspension") ///
	xlabel(-6 "Jul 24" -3 "Oct 24" 0 "Jan 25" 3 "Apr 25" ///
		6 "Jul 25" 9 "Oct 25" 12 "Jan 26") ///
	graphregion(color(white)) plotregion(margin(small)) ///
	legend(order(2 "Region + Month FE (TWFE)" 4 "Region + Country×Month FE") ///
		rows(1) size(small) position(6))
graph export "$figures/twfe_vs_cym_events.png", replace width(1200)

restore

di _newline "Done — subnational event study plots (TWFE) + two-model comparison."
