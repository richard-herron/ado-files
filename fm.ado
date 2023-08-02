*! Date		:	2020-06-26
*! Version	:	0.15
*! Author	:	Richard Herron
*! Email	:	richard_herron@icloud.com

*! takes coefficients from -statsby- and generates Newey-West SEs

/*
2020-06-26 v0.15 added time-series operations and removed FORCE and xtewreg
2019-05-05 v0.14 option to use Toni Whited's xtewreg
2018-03-19 v0.13 option to save coefficients
2018-02-11 v0.12 allow dummies
2017-12-08 v0.11 after preserve, aggressively subset for speed
2017-11-08 v0.10 force option to force irregular time series
2017-09-07 v0.9 allow abbreviation of options
2017-06-30 v0.8 removed marginal effects
2017-06-30 v0.7 marginal effects use sample bhat
2017-06-29 v0.6 logit/probit models return exp(beta*x) marginal effects
2017-06-28 v0.5 marginal effect options (cross-sectional iqr and sd)
2016-12-11 v0.4 unique name for average R2
2016-07-21 v0.3 option to save first-stage results
2016-07-20 v0.2 more flexible, allows arbitrary first-stage regression
2016-07-18 v0.1 first upload to GitHub
*/

program define fm, eclass 
	version 13

	syntax varlist(ts numeric) [if] [in] [ , Estimator(string) Lags(integer 0) Options(string) ]
	marksample touse
	tempname beta VCV
	
	* regress is default estimator
	if "`estimator'" == "" local estimator "regress"

	* add comma prefix to options
	if "`options'" != "" local options ", `options'"

	* get panel variables
	quietly xtset
	local panel `r(panelvar)'
	local time `r(timevar)'

	* parse estimator, y, and X
	tokenize `varlist'
	local y `1'
	macro shift 1
	local X `*'

	preserve

	* estimate first-stage (cross-sectional) coefficients
	if inlist("`estimator'", "regress", "areg") {
		quietly statsby _b N = e(N) r2 = e(r2), by(`time') clear : `estimator' `y' `X' `options'
		rename (_eq2_N _eq2_r2) (N r2)
	}
	else if inlist("`estimator'", "probit", "logit", "logistic") {
		quietly statsby _b N = e(N) r2 = e(r2_p), by(`time') basepop(_n < 10000) clear : `estimator' `y' `X' `options'
		rename (_eq2_N _eq2_r2) (N r2)
	}
	else if inlist("`estimator'", "tobit") {
		quietly statsby _b N = e(N) r2 = e(r2_p), by(`time') basepop(_n < 10000) clear : `estimator' `y' `X' `options'
		rename (_eq3_N _eq3_r2) (N r2)
		capture drop sigma_b_cons
	}
	else {
		display as error "Estimator `estimator' not supported"
		exit 111
	}


	* estimate time-series means and standard errors

	* independent variable SEs first
	quietly tsset `time'
	quietly ds `time' N r2, not
	foreach x of varlist `r(varlist)' {
		if (`lags' > 0) {
			quietly newey `x', lag(`lags')
		}
		else {
			quietly regress `x'
		}
		matrix `beta' = nullmat(`beta'), e(b)
		matrix `VCV' = nullmat(`VCV'), e(V)
	}

	* generate covariance matric from row vector
	matrix `VCV' = diag(`VCV')

	* assign matrix names
	matrix colnames `beta' = `X' _cons
	matrix colnames `VCV' = `X' _cons
	matrix rownames `VCV' = `X' _cons
	* matrix list `beta'
	* matrix list `VCV'

	* generate number of observations and panels
	summarize N, meanonly
	local N = r(sum)
	local T = r(N)
	local df_r = `T' - 1
	local df_m = colsof(`VCV')

	* generate average R-squared
	summarize r2, meanonly
	local r2_avg = r(mean)

	* post results 
	* depname(`y') option requires y to be available
	ereturn post `beta' `VCV', depname("`y'") obs(`N') 
	ereturn scalar df_m = `df_m'
	ereturn scalar df_r = `df_r'
	ereturn scalar T = `T'
	ereturn scalar r2_avg = `r2_avg'

	* return average J-stat for Whited regressions
	if inlist("`estimator'", "xtewreg") {
		ereturn scalar Jstat_avg  = `Jstat_avg'
		ereturn scalar Jstat_p_avg  = `Jstat_p_avg'
	}

	ereturn local cmd "fm"
	if (`lags' > 0) {
		ereturn local vce "Newey-West (1987) standard errors with `lags' lag"
		local title "Fama-Macbeth (1973) regression with Newey-West (1987) standard errors (`lags' lag)"
		ereturn local title `title'
	}
	else {
		ereturn local vce "Fama-Macbeth (1973) standard errors"
		local title "Fama-Macbeth (1973) regression with Fama-Macbeth (1973) standard errors"
		ereturn local title `title'
	}

	* F-test, after posting results
	quietly test `X'
	ereturn scalar F = r(F)
	ereturn scalar p = fprob(e(df_m), e(df_r), e(F))

	* display results
	display as text "`title'"
	display _column(42) as text "First-stage estimator is `estimator'"
	display _column(42) as text "Number of observations"	_column(67) " = " as result %9.0gc e(N)
	display _column(42) as text "Number of cross sections"	_column(67) " = " as result %9.0gc e(T)
	display _column(42) as text "F(" %2.0f e(df_m) ", " %4.0f e(df_r) ")"	_column(67) " = " as result %9.3gc e(F)
	display _column(42) as text "Prob > F"					_column(67) " = " as result %9.3f fprob(e(df_m), e(df_r), e(F))
	display _column(42) as text "Average R-squared"		_column(67) " = " as result %9.3f e(r2_avg)
	if inlist("`estimator'", "xtewreg") {
		display _column(42) as text "Average J-stat"		_column(67) " = " as result %9.3f e(Jstat_avg)
		display _column(42) as text "Average p(J-Stat)"		_column(67) " = " as result %9.3f e(Jstat_p_avg)
	}
	ereturn display

end
