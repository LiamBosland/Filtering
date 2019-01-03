function [var_eps] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim)
% This MATLAB routine closes the loop and computes the variance of the
% residual wavefront var_eps, from system matrices A, G, H, covariance
% matrices Cphi0, sigmae*I, Cw, Kalman-gain matrix K and the wavefront data
% phiSim

% Using the given system matrices, a state-space model can be described
sys_VAR = ss(A,[AH -H],G,[],-1);
% The one
end

