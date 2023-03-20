%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #1
% MAIN Script file
% ECON630 FALL 2022
% Authors: Giuliano Simoncelli, Daniel Schwindt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Call Q1 Script which solves minimization problem using Newton's Method
Q1_PS1;       % Contains functions and code to run particular initial values
heatmap_chart; % Creates Figure 2 in PS1 PDF report
plot3d;       % Creates Figure 1 in PS1 PDF report


%% Call Q2 Script which computes moments from bivariate normal distribution
Q2_PS1;       % Computes moments (i) analytically (ii) Monte Carlo Methods
              % (iii) Gauss-Hermite Quadrature (GHQ) (iv) Tauchen (1986)
