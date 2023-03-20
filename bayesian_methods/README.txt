%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #4
% ECON630 FALL 2022
% Author: Giuliano Simoncelli & Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MASTER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%
PATH: ./Matlab files/
'PS4_Master.m':
Calls the scripts that solve each question


QUESTION 1
%%%%%%%%%%%%%%%%%%%%%%%%%
PATH: ./Matlab files/
1- q1_main.m:
Performs maximum likelihood estimation of linear regression model for three different error specifications:
(a) N(0,1)
(b) t with 1 degree of freedom
(c) Laplace with mu=0 and b=1
Additionally creates relevant charts and tables

Support programs:
1- laplace.m:
Function to draw iid laplacian errors for a given mu and b


QUESTION 2
%%%%%%%%%%%%%%%%%%%%%%%%%
PATH: ./Matlab files/
1- rwmh_gibbs.m: 
Performs random draws from distribution using 
(a) Random Walk Metropolis-Hastings Algorithm (function calls from rwmh.m)
(b) Gibbs Sampling

Support Programs:
1- rwmh.m: 
Function to perform number of draws from distribution using Random Walk Metropolis-Hastings Algorithm

Output:
Charts and CSV tables saved to ../output/
1 - x_moments.csv: select estimated moments for x_1 rv
2 - y_moments.csv: select estimated moments for x_2 rv
3 - cont_gibbs.eps: 2D histogram contour plot for X using Gibbs Sampling
4 - cont_rwmh.eps: 2D histogram control plot for X using RWMH algorithm
5 - hist_gibbs.eps: 3D histogram for X using Gibbs Sampling
6 - hist_rwmh.eps: 3D histogram for X using RWMH
