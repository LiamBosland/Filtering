function [var_eps,VAF] = AOloopSID(G,H,As,Cs,Ks,sigmae,phi)
% This MATLAB routine closes the loop for the Subspace Identification model
% and computes the variance of residual wavefront, using system matrices
% As, Cs, Ks, G, H, covariance matric sigmae^2*I and simulated data phiSim

% Define some shorthards
o = size(G,1);  % The amount of outputs per time-step k (equal to 2*p^2)
% Cphi0 = Cphi(phi);
% Cphi0 = 0.5*(Cphi0 + Cphi0');
% Gamma = Cphi0*G'*inv(G*Cphi0*G' + sigmae^2*eye(o));
Gamma = pinv(G);

% Next, the closed-loop is computed for the provided phiSim
[m2,N] = size(phi);     % The size of phi(k) is (m^2)x1 (real), m = p+1
[u,eps_est,eps_true,phi_est] = deal(zeros(m2,N));
[s,s_est] = deal(zeros(o,N));
x_est = zeros(size(As,1),N);
phi(:,1) = phi(:,1) - mean(phi(:,1));
for k = 2 : N
    e = sigmae^2*eye(o)*randn(o,1);
    s(:,k) = G*phi(:,k) + e;
    
    x_est(:,k+1) = (As-Ks*Cs)*x_est(:,k) + Ks*s(:,k);
    s_est(:,k+1) = Cs*x_est(:,k+1);
    phi_est(:,k+1) = Gamma*s_est(:,k+1);
    u(:,k) = inv(H)*Gamma*s_est(:,k+1);
    eps_est(:,k+1) = Gamma*s_est(:,k+1) - H*u(:,k);
    eps_est(:,k+1) = eps_est(:,k+1) - mean(eps_est(:,k));
    eps_true(:,k) = phi(:,k) - H*u(:,k-1);
end

VAF = vaf(phi,phi_est(:,1:end-1));
var_eps = mean(var((eps_true-eps_est(:,1:end-1))));
