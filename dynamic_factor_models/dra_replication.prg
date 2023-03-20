'--------------------------------------------------------------------------
' ECON630: Problem Set #5
' Part C
' Replicate state space model for yield curve in Diebold, Rudebusch, & Aruoba (2006)
' Author: Daniel Schwindt
' Date: November 14, 2022

'--------------------------------------------
' 0. Preliminaries
%path = @runpath
cd %path

'--------------------------------------------
' 1. Load Data
wfopen "..\data\dra.wf1"

'-------------------------------------------
' 2. Estimate state space model
dra_sspace.ml(maxit=500)
dra_sspace.makestates(t=smooth, n=state_series) *f ' make state series

'-----------------------------------------
' 3. Graph output
graph stateseries.line state_series
stateseries.setelem(1) legend("Level")
stateseries.setelem(2) legend("Slope")
stateseries.setelem(3) legend("Curvature")

stateseries.save(t=postscript, u=in, w=10.5, h=7.5, land, -box) "..\output\stateseries"

close @wf


