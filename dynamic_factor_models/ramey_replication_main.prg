'------------------------------------------------------------------------------------------------
' ECON630: Problem Set #5, Part A MAIN
' Replicating Ramey (2016) Figure 3.1, 3.2A, 3.2B
' Author: Daniel Schwindt
' Date: 11/16/2022
' Date updated: 11/18/2022
'------------------------------------------------------------------------------------------------

'%path = @runpath
'cd %path

' Subprograms for CEE (1999) using Bootstrapping (90% CI)
exec ./ramey_replication_cee_allv_e ' All variables 1965-1995
exec ./ramey_replication_cee_allv_l ' All Variables 1983-2007
exec ./ramey_replication_cee_nocm_e ' W/O Commod Price, 1965-1995
exec ./ramey_replication_cee_nocm_l ' W/W Commod Price, 1983-2007

' Subprograms for Coibion (2012) using Eviews built-in bootstrapping
exec ./ramey_replication_coi

' Subprograms for Jorda (2005) Local Projections technique
exec ./ramey_replication_jorda


