'------------------------------------------------------------------------------------------------
' ECON630: Problem Set #5, Part A (1965 - 1995, No Commodity Price)
' Replicating Ramey (2016) Figure 3.1(1965 - 1995, No Commodity Price)
' Author: Daniel Schwindt
' Date: 11/16/2022
' Date updated:
'------------------------------------------------------------------------------------------------

'-----------------------------------------------------------
' 0. Preliminaries
%path = @runpath
cd %path

wfopen "..\data\ramey_data.wf1"

!n_irf = 48 ' set impulse response length

'---------------------------------------------------------
' Christiano, Eichenbaum, Evans (1999)
'--------------------------------------------------------

'----------------------------------------------------------
' 1. Estimate VARs
'----------------------------------------------------------
' W/O Commodity Price Index, 1965m1-1995m6
smpl 1965m1 1995m6
var var_cee_nocm_early.ls 1 12 lip unemp lcpi ffr lnbr ltr lm1
var_cee_nocm_early.impulse(!n_irf, rmat=imp_nocm_e) @imp ffr
var_cee_nocm_early.makeresids(n=res_cee_nocm_early)
'var_cee_nocm_early.impulse(!n_irf, se=boot, bs=sp, rep=500, cilevels=0.90) @imp ffr

'----------------------------------------------------------------------
' 2. Bootstrapping
'---------------------------------------------------------------------
' Convert residual groups into matrices
matrix mres_allv_e = @convert(res_cee_nocm_early)

' Create coefficient matrices
matrix(7,85) coef_cee_allv_early

' Loop to fill in coefficient matrices
for !i = 1 to 7
	for !j = 1 to 85
		coef_cee_allv_early(!i, !j) = var_cee_nocm_early.c(!i, !j)
'		coef_cee_allv_late(!i, !j) = var_cee_allv_late.c(!i, !j)
	next
next

' Create matrices for pre-pending to matrix of bootstrap time series for VAR estimation
matrix(12,7) init_allv_e
matrix(12,7) init_allv_l
matrix(12,7) init_nocm_e
matrix(12,7) init_nocm_l

' Create matrices for bootstrap bands
for %var lip unemp lcpi  ffr lnbr ltr lm1
	matrix(!n_irf,500) imp_allv_e_{%var}
next

for !j= 1 to 500

	' Construct random samples from matrices of residuals
	matrix b_allv_e = @resample(mres_allv_e)

	' Initialize variable vectors
	%span = "1964m01 1964m12"
	matrix v_lip = @convert(lip, %span)
	matrix v_unemp = @convert(unemp, %span)
	matrix v_lcpi = @convert(lcpi, %span)
	matrix v_ffr = @convert(ffr, %span)
	matrix v_lnbr = @convert(lnbr, %span)
	matrix v_ltr = @convert(ltr, %span)
	matrix v_lm1 = @convert(lm1, %span)
	matrix(1) const = 1

	' Fill init matrix
	!c = 1
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		matplace(init_allv_e, v_{%var}, 1,!c)
	!c = !c + 1
	next

	' Re-order variables since coeffs are in diff order
	vector(12) ord
	ord.fill 12,11,10,9,8,7,6,5,4,3,2,1
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		v_{%var} = v_{%var}.@row(ord)
	next

	' Create Right Hand side vector
	matrix(85,1) rhs_allv
	!c = 0
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		matplace(rhs_allv, v_{%var}, !c*@rows(v_{%var})+1, 1)
	!c = !c + 1
	next
	matplace(rhs_allv, const, 7*@rows(v_lip)+1,1)
	' Create matrix to store predicted values 
	matrix(366, 7) new = 0
	matrix boot = @subextract(b_allv_e, 1, 1,1)
	matrix(1,7) tmp = rhs_allv.@t*coef_cee_allv_early.@t + boot
	matplace(new, tmp, 1,1)

for !i = 2 to 366
' Create new matrices with rolling window
	!c = 1
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		matrix v_{%var}_new = new.@sub(!i-1,!c)
		v_{%var} = @subextract(v_{%var}, 1, 1, 11)
		v_{%var} = @vcat(v_{%var}_new, v_{%var})
	!c = !c + 1
	next

	' Create new rhs vector
	!c = 0
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		matplace(rhs_allv, v_{%var}, !c*@rows(v_{%var})+1, 1)
	!c = !c + 1
	next
	'matplace(rhs_allv, const, 7*@rows(v_lip)+1,1)

	' Compute new right hand side variable
	matrix boot = @subextract(b_allv_e, !i, 1, !i)
	tmp = rhs_allv.@t*coef_cee_allv_early.@t + boot
	matplace(new, tmp, !i, 1)

	next

' Pre-pend prior year data
new = @vcat(init_allv_e, new)
' Convert matrix into time series
series lip_t
series unemp_t
series lcpi_t
series ffr_t
series lnbr_t
series ltr_t
series lm1_t
group gt lip_t unemp_t lcpi_t ffr_t lnbr_t ltr_t lm1_t
sample st 1964m1 1995m6
mtos(new, gt, st)

'Restrict sample to 1965m1 to 1995m6
smpl 1965m1 1995m6
var var_tmp.ls 1 12 lip_t unemp_t lcpi_t ffr_t lnbr_t ltr_t lm1_t
var_tmp.impulse(!n_irf, rmat=imp_t) @imp ffr_t

'Fill bootstrap impulse matrices
	!h = 1
	for %var lip unemp lcpi  ffr lnbr ltr lm1
		matrix imp_t_{%var} = @subextract(imp_t, 1, !h, @rows(imp_t), !h)
		matplace(imp_allv_e_{%var}, imp_t_{%var}, 1,!j)
		!h = !h + 1
	next

delete imp_t st

next

' Rowsort matrices
for %var lip unemp lcpi  ffr lnbr ltr lm1
	imp_allv_e_{%var} = @rowsort(imp_allv_e_{%var}, "a")
next

' Pick top and bottom 50th cols
vector(2) pick
pick.fill 50, 450
!i =1
for %var lip unemp lcpi  ffr lnbr ltr lm1
	matrix bands_imp_allv_e_{%var} = imp_allv_e_{%var}.@col(pick)
	matrix imp_allv_e_{%var}_f = @subextract(imp_nocm_e, 1,!i,@rows(imp_nocm_e), !i)
	imp_allv_e_{%var}_f = @hcat(imp_allv_e_{%var}_f, bands_imp_allv_e_{%var})
	!i = !i + 1
next

' Create graphs
for %var lip unemp lcpi  ffr lnbr ltr lm1
	freeze(gr_{%var}) imp_allv_e_{%var}_f.line
	gr_{%var}.setelem(2) lcolor(orange)
	gr_{%var}.setelem(3) lcolor(orange)
	gr_{%var}.setelem(2) lwidth(1)
	gr_{%var}.setelem(3) lwidth(1)
	gr_{%var}.legend -display
	gr_{%var}.addtext(t,ac) Response of {%var} to FFR
next

graph gr_cee_nocm_e.merge gr_lip gr_unemp gr_lcpi gr_ffr gr_lnbr gr_ltr gr_lm1
gr_cee_nocm_e.align(3, 1, 1.5)
gr_cee_nocm_e.addtext(t,ac) Response to Cholesky One S.D. (d.f. adjusted) Innovations \r CI using Standard percentile bootstrap with 500 bootstrap reps
gr_cee_nocm_e.save(t=postscript, u=in, w=10.5, h=7.5, land, -box) "..\output\gr_cee_nocm_e"
	'matrix b_allv_l = @resample(mres_allv_l)
	'matrix b_nocm_e = @resample(mres_nocm_e)
	'matrix b_nocm_l = @resample(mres_nocm_l)
close @wf


