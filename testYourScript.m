clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                          CAUTION                                            %%
%% This test will only help you spot syntax errors and check for dimension mismatches.         %%
%% Passing the test does not mean your script will run error free with another flight dataset. %%
%% It is YOUR responsibility to design your script to ROBUSTLY handle different scenarios.     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data to test your design
% load('dataTask3_1.mat') %with only input measurement sensor faults
% load('dataTask3_2.mat') %with only angle of attack sensor fault
load('dataTask4.mat') %with both input measurement sensor and angle of attack sensor faults.

%Test your script, remember to change the function name according to your student ID (as on Blackboard)
tic;
[x_est,b_est,Ax_f_instance,Ay_f_instance,Az_f_instance,p_f_instance,q_f_instance,r_f_instance,AoA_f_instance] = integrated_navigation(c_k, d_k, t, dt);
tc=toc;

% Please make sure your code finishes running within 5 minutes.
disp(['Total run time is ',num2str(tc),' seconds'])

if size(x_est,1)==length(t) && size(x_est,2)==12 && ~any(isnan(x_est),'all')
    disp('Dimension of x_est returned by the function is in good order')
else
    error('Dimension of x_est returned by the function is not correct or x_est contains NaN elements')
end

if size(b_est,1)==length(t) && size(b_est,2)==6 && ~any(isnan(b_est),'all')
    disp('Dimension of b_est returned by the function is in good order')
else
    error('Dimension of b_est returned by the function is not correct or b_est contains NaN elements')
end

