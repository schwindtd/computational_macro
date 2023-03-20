ECON630 PS5 README

PART A

File Structure
./ : parent directory
./code/: contains below program files
./data/: contains data files
./output/: contains output files (e.g., charts)

List of Programs
1 - ramey_replication_main.prg: MAIN program calling the following sub-programs:
1a- ramey_replication_cee_allv_e.prg: replicates Figure 3.1 Ramey (2016) IRFs for CEE (1999) using bootstrap CIs, ALL variables & 1965m1-1995m6
1b- ramey_replication_cee_allv_l.prg: same as above but for 1983m1-2007m12
1c- ramey_replication_cee_nocm_e.prg: same as 1a but w/o commodity price variable in VAR
1d- ramey_replication_cee_nocm_l.prg: same as 1b but w/o commodity price variable in VAR
1e- ramey_replication_coi.prg: reproduces Figure 3.2A Ramey (2016) IRFs using Eviews built-in bootstrap CIs
1f- ramey_replication_jorda.prg: reproduces Figure 3.2B Ramey (2016) IRFs using Jorda (2005) local projections

List of Data
1- ramey_data.wf1: Eviews workfile dataset created by importing Monthly tab from Monetarydat.xlsx (sourced from Ramey Website)

List of Output
1- gr_cee_allv_e.eps: CEE (1999) IRFs | all variables | 1965m1-1995m6 | bootstrap 90% CI
2- gr_cee_allv_l.eps: CEE (1999) IRFs | all variables | 1983m1-2007m12 | bootstrap 90% CI
3- gr_cee_nocm_e.eps: CEE (1999) IRFs | w/o commod price | 1965m1-1995m6 | bootstrap 90% CI
4- gr_cee_nocm_l.eps: CEE (1999) IRFs | w/o commod price | 1983m1-2007m12 | bootstrap 90% CI
5- gr_coibion_e.eps: Coibion (2012) IRFs 1965m1-1995m6 bootstrap 90% CI
6- gr_coibion_l.eps: Coibion (2012) IRFs 1983m1-2007m12 bootstrap 90% CI
7- gr_jorda.eps: Jorda (2005) Local Projections IRFs 1969m3- 1996m12 90% CI

PART B

File Structure
./ : parent directory
./code/: contains below program files
./data/: contains data files
./output/: contains output files (e.g., charts)

List of Programs
1- PartB.prg: Replicate Illustration 2.1 from Del Negro & Schorfheide (2011)

List of Data
1- PartB_workfile.wf1: Eviews workfile created by fetching data from the FRED database

List of Output
1- partb_percentile.png: Replication of Fig.2 in Del Negro & Schorfheide (2011), with the 5th and 95th percentiles bands
2- partb_ci.png: Replication of Fig.2 in Del Negro & Schorfheide (2011), with the 90% confidence interval bands

PART C
File Structure
./ : parent directory
./code/: contains below program files
./data/: contains data files
./output/: contains output files (e.g., charts)

List of Programs
1 - dra_replication.prg: Replicates state space model for yield curve in Diebold, Rudebusch & Aruoba (2006)

List of Data
1- dra.wf1: Eviews workfile dataset created by importing DRA Data.txt file from Aruoba website

List of Output
1- stateseries.emf
2- Estimated Coefficients: obtained from dra_sspace object in workfile



