function [var_eps] = AOloopSID(G,H,As,Cs,Ks,sigmae,phiSim)
% This MATLAB routine closes the loop for the Subspace Identification model
% and computes the variance of residual wavefront, using system matrices
% As, Cs, Ks, G, H, covariance matric sigmae^2*I and simulated data phiSim

% Define some shorthards
B = [As*H -H];
o = size(G,1);  % The amount of outputs per time-step k (equal to 2*p^2)

% Next, the closed-loop is computed for the provided phiSim
[m2,N] = size(phiSim);     % The size of phi(k) is (m^2)x1 (real), m = p+1
[u,eps_est,eps_true,phi_DM] = deal(zeros(m2,N));
s = [(G*phiSim(:,1) + (sigmae*eye(o)*randn(o,1))) zeros(o-1,N)];
eps_est(:,2) = B*[zeros(m2,1) ; u(:,2)] + Ks*(s(:,2) - Cs*Gamma*s_est(:,2));
for k = 2 : N
    e = sigmae^2*eye(o)*randn(o,1);
    phi_DM(:,k) = H*u(:,k-1);
    eps_true(:,k) = phiSim(:,k) - phi_DM(:,k);
    eps_true(:,k) = eps_true(:,k) - mean(eps_true(:,k));
    s(:,k) = G*eps_true(:,k) + e;
    u(:,k) = inv(H)*(As*eps_est(:,k) + As*H*u(:,k-1) + K*(s(:,k) - Cs*Gamma*s_est(:,k)));
    eps_est(:,k+1) = As*eps_est(:,k) + B*[u(:,k-1) ; u(:,k)] + Ks*(s(:,k) - Cs*Gamma*s_est(:,k));
end

