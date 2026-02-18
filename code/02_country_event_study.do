/*
02_country_event_study.do
Country-level event study coefplot (re-centred: pre-period mean = 0).
Produces: event_study_armed.png
Requires: $output/acled_conflict_panel.dta
*/

clear all
set more off

// SET THIS TO YOUR PROJECT ROOT
global root "/Users/lee/Library/CloudStorage/Dropbox-CGDEducation/Lee Crawfurd/blogs/2026-02-usaid-conflict"

global input   "$root/data/input"
global output  "$root/data"
global figures "$root/figures"

use "$output/acled_conflict_panel.dta", clear
keep if event_time >= -6 & event_time <= 12

// Event-time interactions (omit t = -1)
forvalues t = -6/12 {
	if `t' == -1 continue
	local tlab = cond(`t' < 0, "m" + string(abs(`t')), "p" + string(`t'))
	gen et_`tlab' = (event_time == `t') * usaid_pc_std
}

local etvars et_m6 et_m5 et_m4 et_m3 et_m2 et_p0 et_p1 et_p2 et_p3 et_p4 et_p5 et_p6 et_p7 et_p8 et_p9 et_p10 et_p11 et_p12

// ============================================================
// Program: re-centred coefplot via twoway
// ============================================================
capture program drop recentred_plot
program define recentred_plot
	syntax, depvar(string) etvars(string) color(string) ///
		title(string) filename(string)

	reghdfe `depvar' `etvars', a(cid ym) cl(cid)

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
		subtitle("Country FE + Month FE, SEs clustered at country") ///
		ytitle("Coefficient (IHS, re-centred)") ///
		xtitle("Months relative to USAID suspension") ///
		xlabel(-6 "Jul 24" -3 "Oct 24" 0 "Jan 25" 3 "Apr 25" ///
			6 "Jul 25" 9 "Oct 25" 12 "Jan 26") ///
		graphregion(color(white)) plotregion(margin(small)) ///
		legend(off) ///
		note("N = `N', Countries = `Ncl'. Filled = p<0.05. Pre-period re-centred to zero.")
	graph export "$figures/`filename'", replace width(1200)

	restore
end

// ============================================================
// Generate plots
// ============================================================
recentred_plot, depvar(ihs_armed) etvars(`etvars') color(cranberry) ///
	title("Armed Conflict — Country-Level Event Study") ///
	filename(event_study_armed.png)

di _newline "Done — event_study_armed.png"
