function [vf,grad_vf] = compute_volume_fraction()
%
% This function computes the volume fraction and its sensitivities
% based on the last geometry projection
%
global FE OPT

% compute the volume fraction
    v_e = sum(OPT.dv*1.875); % element volumes
    V = FE.n_elem; % total volume
    vf = v_e/V;

% compute the design sensitivity
    grad_vf =OPT.H'*ones(FE.n_elem,1)/V;
    
% output
    OPT.volume_fraction = vf;
    OPT.grad_volume_fraction = grad_vf;