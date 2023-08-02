* sample usage: bootstrap, cluster(firm) idcluster(firm2) : xtpoissoncf `v' $FF01 $CONTROLS i.year, endog(gov41 gov41_antisd) exog(gov41_iv gov41_iv_antisd) panel(firm2 year) resid(resid)

program define xtpoissoncf, eclass
	version 13

	syntax varlist(fv) [if] [in], endog(varlist) exog(varlist) panel(varlist) resid(string)
	marksample touse

	* parse varlists
	tokenize `varlist'
    local y `1'
    macro shift
    local X `*'
    tokenize `panel'
    local i `1'

	* generate first-stage residuals
	xtset `panel'
	local j = 1
	foreach v of local endog {
		xtreg `v' `exog' `X' if `touse', fe
		capture drop `resid'`j'
		predict `resid'`j', e
		local ++j 
	}

	* estimate fixed effects poisson
	xtpoisson `y' `endog' `resid'* `X' if `touse', fe
	tempname b VCV
	matrix `b' = e(b)
	matrix `VCV' = e(V)
	local N = e(N)
	replace `touse' = e(sample)

	* weak IV tests from xtivreg2
	xi: xtivreg2 `y' `X' (`endog' = `exog') if `touse', fe robust cluster(`i')
	local cdf = e(cdf)
	local rkf = e(rkf)

	* return for estout/esttab
	ereturn post `b' `VCV', obs(`N') esample(`touse') depname("`y'")
	ereturn scalar cdf = `cdf'
	ereturn scalar rkf = `rkf'
	ereturn display
end


