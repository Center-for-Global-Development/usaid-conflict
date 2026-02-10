/*
04_did_coefficients.do
All difference-in-differences regression coefficients.
Main spec: TWFE (admin1 + month FE).
Robustness: admin1 + country×month FE, and + admin1-specific trends.
*/

clear all
set more off

// SET THIS TO YOUR PROJECT ROOT
global root "/Users/lee/Library/CloudStorage/Dropbox-CGDEducation/Lee Crawfurd/blogs/2026-02-usaid-conflict"

global input   "$root/data/input"
global output  "$root/data"
global figures "$root/figures"

di _newline "============================================"
di "ALL DiD COEFFICIENTS"
di "============================================"

// ============================================================
// 1. COUNTRY-LEVEL (country + month FE)
// ============================================================
use "$output/acled_conflict_panel.dta", clear
keep if event_time >= -6 & event_time <= 12
gen post_treat = post * usaid_pc_std

gen ihs_fat_armed = asinh(fat_armed)
gen ihs_fat_all = asinh(fat_all)

di _newline "--- COUNTRY-LEVEL (country + month FE, clustered at country) ---"

reghdfe ihs_armed post_treat, a(cid ym) cl(cid)
di "IHS armed events:      coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_all post_treat, a(cid ym) cl(cid)
di "IHS all events:        coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_fat_armed post_treat, a(cid ym) cl(cid)
di "IHS armed fatalities:  coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_fat_all post_treat, a(cid ym) cl(cid)
di "IHS all fatalities:    coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

// ============================================================
// 2. SUBNATIONAL — TWFE (admin1 + month FE) [MAIN SPEC]
// ============================================================
use "$input/subnational_panel.dta", clear
gen ihs_fat_armed = asinh(fat_armed)
gen ihs_fat_all = asinh(fat_all)

keep if event_time >= -6 & event_time <= 12
encode region_id, gen(rid)
egen cym = group(country ym)
gen post_treat = post * usaid_std

di _newline "--- SUBNATIONAL TWFE (admin1 + month FE, clustered at admin1) [MAIN] ---"

reghdfe ihs_events_armed post_treat, a(rid ym) cl(rid)
di "IHS armed events:      coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_events_all post_treat, a(rid ym) cl(rid)
di "IHS all events:        coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_fat_armed post_treat, a(rid ym) cl(rid)
di "IHS armed fatalities:  coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_fat_all post_treat, a(rid ym) cl(rid)
di "IHS all fatalities:    coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

// ============================================================
// 3. ROBUSTNESS: admin1 + country×month FE
// ============================================================
di _newline "--- ROBUSTNESS: admin1 + country×month FE ---"

reghdfe ihs_events_armed post_treat, a(rid cym) cl(rid)
di "IHS armed events:      coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

reghdfe ihs_fat_armed post_treat, a(rid cym) cl(rid)
di "IHS armed fatalities:  coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat])) "  N = " e(N) "  cl = " e(N_clust)

// ============================================================
// 4. ROBUSTNESS: admin1 + country×month FE + admin1 linear trends
// ============================================================
di _newline "--- ROBUSTNESS: + admin1 linear trends ---"

reghdfe ihs_events_armed post_treat, a(rid cym rid#c.event_time) cl(rid)
di "IHS armed events:      coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat]))

reghdfe ihs_fat_armed post_treat, a(rid cym rid#c.event_time) cl(rid)
di "IHS armed fatalities:  coef = " %7.4f _b[post_treat] "  se = " %7.4f _se[post_treat] "  p = " %5.3f 2*ttail(e(df_r), abs(_b[post_treat]/_se[post_treat]))

di _newline "Done."
