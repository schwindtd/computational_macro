'------------------------------------------------------------------------------------------------
' ECON630: Problem Set #5, Part A
' Replicating Ramey (2016) Figure 3.2B
' Author: Daniel Schwindt
' Date: 11/15/2022
' Date updated:
'------------------------------------------------------------------------------------------------

'-----------------------------------------------------------
' 0. Preliminaries
%path = @runpath
cd %path

wfopen "..\data\ramey_data.wf1"

!n = 48

'--------------------------------------------------------------------------------
' Jorda (2005) Local Projection with Romer & Romer Shocks
'--------------------------------------------------------------------------------
smpl 1969m3 1996m12
' Estimate VAR
var var_jorda.ls 1 12 lip unemp lcpi rrshock
var_jorda.makeresids(n=res_jorda)

' Loop to perform Local Projection per Ramey (2016)
matrix(!n+1,4) jorda_coeffs
matrix(!n+1,4) jorda_ses
for !j = 0 to !n
equation jorda_lip.ls(cov=hac, covbw=!j+1) lip(!j) rrshock rrshock(-1) rrshock(-2) ffr(-1) ffr(-2) lip lip(-1) lip(-2) unemp unemp(-1) unemp(-2) lcpi lcpi(-1) lcpi(-2) lpcom lpcom(-1) lpcom(-2)  
jorda_coeffs(!j+1,1) = jorda_lip.@coef(1)
jorda_ses(!j+1,1) = jorda_lip.@stderrs(1)

equation jorda_unemp.ls(cov=hac, covbw=!j+1) unemp(!j) rrshock rrshock(-1) rrshock(-2) ffr(-1) ffr(-2) lip lip(-1) lip(-2) unemp unemp(-1) unemp(-2) lcpi lcpi(-1) lcpi(-2) lpcom lpcom(-1) lpcom(-2)  
jorda_coeffs(!j+1,2) = jorda_unemp.@coef(1)
jorda_ses(!j+1,2) = jorda_unemp.@stderrs(1)

equation jorda_lcpi.ls(cov=hac, covbw=!j+1) lcpi(!j) rrshock rrshock(-1) rrshock(-2) ffr(-1) ffr(-2) lip lip(-1) lip(-2) unemp unemp(-1) unemp(-2) lcpi lcpi(-1) lcpi(-2) lpcom lpcom(-1) lpcom(-2)  
jorda_coeffs(!j+1,3) = jorda_lcpi.@coef(1)
jorda_ses(!j+1,3) = jorda_lcpi.@stderrs(1)

equation jorda_ffr.ls(cov=hac, covbw=!j+1) ffr(!j) rrshock rrshock(-1) rrshock(-2) ffr(-1) ffr(-2) lip lip(-1) lip(-2) unemp unemp(-1) unemp(-2) lcpi lcpi(-1) lcpi(-2) lpcom lpcom(-1) lpcom(-2)  
jorda_coeffs(!j+1,4) = jorda_ffr.@coef(1)
jorda_ses(!j+1,4) = jorda_ffr.@stderrs(1)
next

' Create Error Bands
matrix(!n+1,4) errs_low = jorda_coeffs - 1.68*jorda_ses
matrix(!n+1,4) errs_high = jorda_coeffs + 1.68*jorda_ses

' Create Matrices for plotting
!i=1
for %var lip unemp lcpi ffr
	matrix(!n+1,3) {%var}_imp
	vector {%var}_errs_low = errs_low.@col(!i)
	vector {%var}_errs_high = errs_high.@col(!i)
	vector {%var}_coef = jorda_coeffs.@col(!i)
	colplace({%var}_imp, {%var}_errs_low, 1)
	colplace({%var}_imp, {%var}_coef,2)
	colplace({%var}_imp, {%var}_errs_high,3)
	!i= !i+1
next

'---------------------------------------------------------
' Chart output
'---------------------------------------------------------
for %var lip unemp lcpi ffr
	freeze(gr_{%var}) {%var}_imp.mixed band(1,3) line(2)
	gr_{%var}.legend -display
	gr_{%var}.addtext(t,ac) {%var}
next

graph gr_full.merge gr_ffr gr_lip gr_unemp gr_lcpi
gr_full.save(t=postscript, u=in, w=10.5, h=7.5, land, -box) "..\output\gr_jorda"

close @wf


