'ECON630-PS5
'Part B: sign restrictions
'Authors: Daniel Schwindt, Luying Yang
'Last edit: 11/22/2022

'---------------Create workfile, Import data-----------------------------------

'wfcreate q 1959Q1 2006Q3
'fetch(link, m) GDPC96 CNP16OV GDPDEF FEDFUNDS  GDP M2SL
%path = @runpath
cd %path

wfopen "..\data\PartB_workfile.wf1"

'---------------Transform  Data--------------------------------------------------
'detrend output per capita
genr outputpc= log(gdpc96/cnp16ov)
equation eq1.ls outputpc c @trend
genr output= resid*100

'inverse velocity
equation eq2.ls log(m2sl/gdp) c @trend
genr  iv = resid*100
genr rmb = iv+output

'inflation
genr infl =400*( log(gdpdef) - log(gdpdef(-1)))

'P matrix
smpl 1965Q1 2005Q1
var mainvar.ls(noconstant) 1 4 output infl fedfunds rmb
matrix cholesky = @cholesky(mainvar.@residcov)
mainvar.impulse(40, se=a,rmat=chol_irf) @imp output infl fedfunds rmb

'---------------Cholesky Responses----------------------------------------------------
!draw = 1000
!period = 40

!k = 0
for %var output infl ffr rmb
	matrix(!period, 4) r_{%var}
	matrix(!period,!draw)  irf_{%var}
	for !i = 1 to 4
     		colplace(r_{%var},chol_irf.@col(!i+!k*4),!i)
	next
	!k = !k+1
next

'---------------Draw Candidates----------------------------------------------------
!val =1
!iter = 1
while !val<=!draw and !iter <50000
''randomly draw alpha
scalar a1 = @nrnd
scalar a2 = @nrnd
scalar a3 = @nrnd
scalar a4 = @nrnd
!norm = (a1^2+a2^2+a3^2+a4^2)^0.5
vector(4) aalpha
aalpha.fill a1/!norm, a2/!norm,  a3/!norm, a4/!norm
if aalpha(3)<0 then
 aalpha = -aalpha
endif
''generate P1
vector cand = cholesky*aalpha
''get irf
for %var output infl ffr rmb
	vector vecirf_{%var} = r_{%var}*cand
next
''test the sign
!sign = (vecirf_infl(1)<0) and  (vecirf_rmb(1)<0) and  (vecirf_ffr(1)>0) and (vecirf_infl(2)<0) and  (vecirf_rmb(2)<0) and  (vecirf_ffr(2)>0) 
if !sign = 1 then
	for %var output infl ffr rmb
		colplace(irf_{%var},vecirf_{%var},!val)
	next
	!val= !val+1
endif
!iter = !iter+1
wend

'---------------Graph the IRF and 90% bands-------------------------------------
for %var output infl ffr rmb
  irf_{%var}=@transpose(irf_{%var})
next

'Get mean irf and bands
for %var output infl ffr rmb
 '90% percentile
  matrix(!period,3) graph_{%var}
  colplace(graph_{%var}, @cquantile(irf_{%var}, 0.95),1)
  colplace(graph_{%var}, @cmean(irf_{%var}),2)
  colplace(graph_{%var}, @cquantile(irf_{%var}, 0.05),3)
  '90% ci
  matrix(!period,3) graph_{%var}_ci
  colplace(graph_{%var}_ci, @cmean(irf_{%var})+1.645*@cstdev(irf_{%var}),1)
 colplace(graph_{%var}_ci, @cmean(irf_{%var}),2)
  colplace(graph_{%var}_ci, @cmean(irf_{%var})-1.645*@cstdev(irf_{%var}),3)
next

'Create graphs
svector(4) names
names(1) = "Output"
names(2) = "Inflation"
names(3) = "Federal Funds Rate"
names(4) = "Real Money Balance"

''IRFs with 5th & 95th percentile bands
!i=1
for %var output infl ffr rmb
     'delete gr_{%var}
	'delete  gr_{%var}_ci
	freeze(gr_{%var}) graph_{%var}.line
	gr_{%var}.setelem(1)  lcolor(red) lwidth(0.8)
	gr_{%var}.setelem(2)  lcolor(black)
	gr_{%var}.setelem(3)  lcolor(red) lwidth(0.8) 
	gr_{%var}.legend -display
	%name = names(!i)
	gr_{%var}.addtext(t,ac,font(20))  {%name}
	!i = !i+1
next
gr_output.axis(left) range(-1,0.6)
gr_infl.axis(left) range(-1,0.6) 
gr_ffr.axis(left) range(-1,1) 
gr_rmb.axis(left) range(-2,1)

graph percentile.merge gr_output gr_infl gr_ffr gr_rmb
percentile.addtext(t,ac,font(28,+b)) "Response  to a Monetary Policy Shock"
percentile.addtext(b,bc,font(20)) "Black lines depict the mean IRF, and red lines correspond to the 5th and 95th percentiles."
show percentile

''IRFs with 90% CI
!i = 1
for %var output infl ffr rmb
	'delete  gr_{%var}_ci
	freeze(gr_{%var}_ci) graph_{%var}_ci.line
	gr_{%var}_ci.setelem(1)  lcolor(red) lwidth(0.8)
	gr_{%var}_ci.setelem(2)  lcolor(black)
	gr_{%var}_ci.setelem(3)  lcolor(red) lwidth(0.8) 
	gr_{%var}_ci.legend -display
	%name = names(!i)
	gr_{%var}_ci.addtext(t,ac,font(20))  {%name}
	!i = !i+1
next

gr_output_ci.axis(left) range(-1,0.6)
gr_infl_ci.axis(left) range(-1,0.6) 
gr_ffr_ci.axis(left) range(-1,1) 
gr_rmb_ci.axis(left) range(-2,1)

graph ci.merge gr_output_ci gr_infl_ci gr_ffr_ci gr_rmb_ci
ci.addtext(t,ac,font(28,+b)) "Response  to a Monetary Policy Shock"
ci.addtext(b,bc,font(20)) "Black lines depict the mean IRF, and red lines correspond to the 90% confidence interval."
show ci


percentile.save(t=png, c, box, port, w=7.5, h = 5.4, u=in, d=480, trans) "..\output\partb_percentile.png"
ci.save(t=png, c, box, port, w=7.5, h = 5.4, u=in, d=480, trans) "..\output\partb_ci.png"

