%% Stochastic Growth Model
% Value Function Iteration Methods
% 
% Daniel Schwindt
% 9/23/2022
% Runs all methods detailed in PS2 for ECON630 - Fall 2022

%% Housekeeping
clear all
close all
clc
%% Create Runtime and Iteration matrices
runtime = nan(1,10);
iters   = nan(1,10);
errs    = nan(2,10);

%% a. Naive VFI
disp('Naive VFI')
Q2_a
runtime(1) = naiveVFItime;
iters(1)   = iter;
errs(:,1)  = [max_err; max_err_pf];
%% b. Normalized Value Function with Steady State as initial guess
disp('Normalized VFI with SS Guess')
Q2_b
runtime(2) = time;
iters(2) = iter;
errs(:,2) = [max_err; max_err_pf];

%% c. Multi-grid n1=50, n2=150, n3=300
disp('Multi-grid method: n1=50, n2=150, n3=300')
Q2_c
runtime(3:6) = [time0, time1, time2, time3];
iters(3:6) = [iter, iter5+iter6, iter7+iter8, iter9+iter10];
errs(:,3:6) = [max_err0, max_err1, max_err2, max_err3; 
               max_err_pf1, max_err_pf1, max_err_pf2, max_err_pf3];

%% d. Multi-grid n1=50 but with monotonicity
disp('Multi-grid n1=50 with monotonicity')
Q2_d
runtime(7) = time;
iters(7) = iter;
errs(:,7) = [max_err; max_err_pf];
%% e. Same as b. but with binary search and monotonicity
disp('Binary Search and Monotonicity')
Q2_e
runtime(8) = time;
iters(8) = iter;
errs(:,8) = [max_err; max_err_pf];
%% g. Endogenous Grid Method
disp('Endogenous Grid Method')
Q2_g
runtime(10) = time;
iters(10) = iter;
errs(:,10) = [max_err; max_err_pf];
%% f. Same as b. but using FOC
disp('FOC')
Q2_f
runtime(9) = time; % takes about an hour to run
iters(9) = iter;
errs(:,9) = [max_err; max_err_pf];
%% Output
header = {'Naive VFI', 'Normalized SS', 'Multi-Grid: 50, 150, 300', ...
          'Multi-Grid: 50', 'Multi-Grid: 150', 'Multi-Grid: 300', ...
          'Multi-Grid: Monoton', 'Binary Search', 'FOC', 'EGM'};
runtime = [header; num2cell(runtime)];
iters = [header; num2cell(iters)];
errs = [header; num2cell(errs)];
writecell(runtime, './output/runtime.csv');
writecell(iters, './output/iters.csv');
writecell(errs, './output/errs.csv');