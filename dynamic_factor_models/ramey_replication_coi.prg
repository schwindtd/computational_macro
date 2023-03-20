'------------------------------------------------------------------------------------------------
' ECON630: Problem Set #5, Part A
' Replicating Ramey (2016) Figures 3.1, 3.2A
' Author: Daniel Schwindt
' Date: 11/15/2022
' Date updated:
'------------------------------------------------------------------------------------------------

'-----------------------------------------------------------
' 0. Preliminaries
%path = @runpath
cd %path

wfopen "..\data\ramey_data.wf1"

'----------------------------------------------------------------------
' Coibion (2012) with Romer & Romer Shocks
'----------------------------------------------------------------------
'----------------------------------------------------------
' 1. Estimate VARs
'----------------------------------------------------------
' Full Coibion variable set, 1965m1-1995m6
smpl 1965m1 1995m6
var var_coi_early.ls 1 12 lip unemp lcpi cumrrshock
'var_coi_early.impulse(50) @imp ffr
'var_coi_early.makeresids(n=res_coi_early)
freeze(gr_e) var_coi_early.impulse(48, se=boot, bs=sp, rep=500, cilevels=0.90) @imp cumrrshock

' Full Coibion variable set, 1983m1-2007m12
smpl 1983m1 2007m12
var var_coi_late.ls 1 12 lip unemp lcpi cumrrshock
'var_coi_late.impulse(50) @imp ffr
'var_coi_late.makeresids(n=res_coi_late)
freeze(gr_l) var_coi_late.impulse(48, se=boot, bs=sp, rep=500, cilevels=0.90) @imp cumrrshock

'----------------------------------------------------------
' 2. Output Charts
'----------------------------------------------------------
gr_e.save(t=postscript, u=in, w=10.5, h=7.5, land, -box) "..\output\gr_coibion_e"
gr_l.save(t=postscript, u=in, w=10.5, h=7.5, land, -box) "..\output\gr_coibion_l"

close @wf


