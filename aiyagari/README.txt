----------------------------------------------------
% Problem Set #6
% ECON630 FALL 2022
% Author: Tom Boesche & Daniel Schwindt
----------------------------------------------------

File Structure
./ : parent directory
./code/: contains below program files
./output/: contains output files (e.g., charts)

------------------------------------------------------
QUESTION 1

List of Programs
1- q1.m: MAIN Q1 program computing SRCE for XI=0 and XI=0.5, calls following sub-programs
1a- srce.m: matlab function that computes Stochastic Recursive Competitive Equilibrium
1a1- pack_mgm_input.m: converts named parameters into one structure to pass to Value Function Iteration sub-program
1a2- mgm_input.m: computes value and policy functions under CRRA utility using multi-grid method
1a3- statdist.m: iterates over (asset, efficiency units) distribution to obtain stationary distribution
1a4- tauchen.m: uses Tauchen (1986) method to discretize AR(1) process

2- asset_demand_supply.m: loops over grid of interest rates to compute asset supply and demand
2a- pack_mgm_input.m: converts named parameters into one structure to pass to Value Function Iteration sub-program
2b- mgm_input.m: computes value and policy functions under CRRA utility using multi-grid meth
2c- tauchen.m: uses Tauchen (1986) method to discretize AR(1) process

List of Outputs
1- Lorenz Curve charts
2- Marginal Distribution of Assets charts
3- Cumulative Distribution of Assets charts
4- Household Decision Rules charts
5- SRCE definitions for XI=0 and XI=0.5

------------------------------------------------------
QUESTION 2

List of Programs
1- q2.m: MAIN Q2 program computing transition path between two SRCEs (Z=1 and Z=0.9)
1a- srce.m: see above description
1a1- pack_mgm_input.m: converts named parameters into one structure to pass to Value Function Iteration sub-program
1a2- mgm_input.m: computes value and policy functions under CRRA utility using multi-grid method
1a3- statdist.m: iterates over (asset, efficiency units) distribution to obtain stationary distribution
1a4- tauchen.m: uses Tauchen (1986) method to discretize AR(1) process

List of Outputs
1- Transition Path charts (r, K, w)