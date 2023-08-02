clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath 'FE_routines'
addpath 'functions'
addpath 'mesh_utilities'
addpath 'optimization'
addpath 'utilities'
addpath 'plotting'
addpath 'Cubic_Spline'
% This function performs a finite difference check of the sensitivities of
% the CONSTRAINT function with respect to the bar design variables.
% It is currently setup for one constraint, but it can be easily modified
% for other/more constraints.
global TPMS FE 
TPMS = 'Primitive';
FE.dim = 2;
tic
init_SGs();
no_FD = 30;
min_t = 0.01;
max_t = 0.05;
dv_0 = min_t + rand(no_FD,1).*(max_t-min_t);
% ===============================
% FINITE DIFFERENCE SENSITIVITIES
% ===============================
n_dv = no_FD;
grad_theta_i = zeros(n_dv, 1);
fd_step = 1e-8;
error = zeros(n_dv,1);
rel_error = zeros(n_dv,1);
grad_theta_FD = zeros(n_dv,1);
grad_theta_Ana = zeros(n_dv,1);
% Finite differences
disp('Computing finite difference sensitivities...');
% Do this for all design variables or only a few
for i = 1:n_dv
    dv_i = dv_0( i ) + fd_step;
    [theta_0,grad_theta_0] = compute_P_norm_SGs(dv_0(i));
    grad_theta_Ana(i,1)  = grad_theta_0;
    [theta_i,~] = compute_P_norm_SGs(dv_i);
    grad_theta_FD(i,1)  = (theta_i - theta_0)/fd_step;
    error(i,1) = grad_theta_Ana(i,1) - grad_theta_FD(i,1);
    rel_error(i,1) = error(i,1)./abs(theta_0);
end
[max_error, ind_max_error] = max(error);
[max_rel_error, ind_max_rel_error] = max(rel_error);
%--------------------
figure();    clf;    hold on
plot(grad_theta_FD,'o','LineStyle','-')
plot(grad_theta_Ana,'.','LineStyle','-')
legend('fd','analytical')
title('Constraint function','Interpreter','latex',  'FontSize', 14, 'FontWeight','bold')
legend('FD','Analytical',  'FontSize', 10, 'FontWeight','bold')
xlabel('Design variable: v', 'FontSize', 12, 'FontWeight','bold')
ylabel('dz/dv',  'FontSize', 12, 'FontWeight','bold')

