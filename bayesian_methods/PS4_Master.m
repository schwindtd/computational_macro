%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Set #4
% MAIN Script file
% ECON630 FALL 2022
% Authors: Daniel Schwindt, Giuliano Simoncelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we investigate the properties of maximum likelihood estimates under
% misspecification. We compute MLE estimates under three different error
% specifications:
% (a) N(0,1)
% (b) t with 1 degree of freedom
% (c) Laplace with mu=0 and b=1

% Script to get all the relevant MLE estimates and charts
%  (Takes about 50 seconds to run)
Q1_MLE
% Output are then stored in '../output/'

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Running this script file will produce random draws from a bivariate
% normal using:
% (a) Random-walk Metropolis-Hastings Algorithm
% (b) Gibbs Sampling
rwmh_gibbs
% Output (charts and CSV tables) are then stored in '../output/'



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
