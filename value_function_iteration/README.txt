%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #2
% ECON630 FALL 2022
% Author: Tom Boesche & Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Compute value function and policy functions, record runtime and iterations and errors vs. true values
1 - Q2_a.m: Naive VFI
2 - Q2_b.m: Normalized Value Function with SS as initial guess
3 - Q2_c.m: Multi-grid n1=50, n2=150, n3=300
    3a - Linear Interpolation: linear_interp.m
    3b - Stochastic VFI function: stoch_gr_vfi.m
4 - Multi-grid n1=50 with monotonicity: Q2_d.m
    4a - Linear Interpolation: linear_interp.m
5 - Q2_e.m: Same as (2) but with monotonicity and binary search: 
6 - Q2_f.m: FOC Method
    6a - Numerical Derivatives: numderiv.m
    6b - First Order Condition: FOC.m
7 - Q2_g.m: Endogenous Grid Method

Each of the above programs 1-7 make use of another program which creates the vectors, matrices, and other variables needed to run the iterations.

A - Environment set-up: stoch_growth_setup.m

In addition, another function file contains the code for the Tauchen (1986) function used in stoch_growth_setup.

B - Tauchen (1986): tauchen.m: written by Jesus Fernandez-Villaverde

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1 - Q3_script.m: 
Computes policy function for new steady state (A=2) using Binary Search/Monotonicity. Then computes the transition path for capital if starting at the prior steady state (A=1). Plots the results.
2 - Q3_script_stoch.m:
Similar to above, but with the TFP shocks not shut off. Used as an additional exercise.