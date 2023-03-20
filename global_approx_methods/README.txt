%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #3
% ECON630 FALL 2022
% Author: Manuel Molina & Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PART A: Code to compute correlation tables across variables in Cooley & Hansen (1995) RBC Model
%%%%%%%%%%%%%%%%%%%%%%%%%
PATH: ./code/dynare/
1- run_linear_ch.m: solves linear, first order RBC model and creates the correlation table
1a- linear_ch.mod: dynare code called by 1- which contains the model and solution steps
2- run_linear_ch_second_order: solves linear, second order RBC model and creates the correlation table
2a- linear_ch_second_order.mod: dynare code called by 2- which contains the model and solution steps
3- run_log_linear_ch.m: solves log-linear, first order RBC model, creates the correlation table and plots the irfs
3a- log_linear_ch.mod: dynare code called by 3- which contains the model and solution steps
4- run_log_linear_second_order.m: solves log-linear, second order RBC model and creates the correlation table
4a- log_linear_ch_second_order.mod: dynare code called by 3- which contains the model and solution steps
5- run_log_linear_multiple_simulations.m: solves log-linear and creates the correlation table calculating the mean
	and standard deviation of the simulations
5a-log_linear_multiple_simulations.mod solves the log-linear, first order RBC model multiple times

Support programs:
1- create_correlation_table.m: called by run_xxxx.m program files in order to create formatted correlation results tables
2- script_solve.m: solves for the steady state of the RBC Cooley & Hansen (1995) model
2a- solve_ss.m: function which defines system of steady state equations, called by script_solve.m
3-  plot_irfs.m called by run_log_linear_ch.m to plot the irfs for the first order log linear model
4- correlation_tables_multiple_simulations.m CALLED BY run_log_linear_multiple_simulations.m to create formatted correlation results tables
PART B
%%%%%%%%%%%%%%%%%%%%%%%%%
PATH: ./code/
1- script_stoch_growth_global.m: main file to compute and plot decision rules for local and global approximations and to plot IRFs of Cooley & Prescott (1995) model with sigma = 20 and sigma_eps=0.035

Support Programs: (Called within script_stoch_growth_global.m)
1- cooley_prescott_linear.mod: dynare file used to compute local approximation of Cooley & Prescott RBC model
2- parameters.m: sets parameter values for script_stoch_growth_global.m
2a- autoreg.m: defines autoregression function for use in setting lower and upper bounds for Z
2b- tauchen.m: function to compute transition probabilities of markov chain approximating AR(1)
3- irf_local.m: function for computing IRFs using local approximation data
4- irf_global.m: function for computing IRFs using global approximation data
