function [var_eps,VAF] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim)
% This MATLAB routine closes the loop of the Vector Auto-Regressive model
% and computes the variance of the residual wavefront var_eps, from system
% matrices A, G, H, covariance matrices Cphi0, sigmae*I, Cw, Kalman-gain
% matrix K and the wavefront data phiSim

% Define some shorthands
B = [A*H -H];
o = size(G,1);  % The amount of outputs per time-step k (equal to 2*p^2)
%y(k) = K*s(k) + (A-K*G)*eps_est(k);

% Next, the closed-loop is computed for the provided phiSim
[m2,N] = size(phiSim);     % The size of phi(k) is (m^2)x1 (real), m = p+1
[y,u,eps_est,eps_true,phi_DM,phi_est] = deal(zeros(m2,N));
s = [(G*eps_est(:,1) + (sigmae*eye(o)*randn(o,1))) zeros(o,N-1)]; % The size of s(k) is 2(p^2)x1 (real)
for k = 2 : N
    e = (sigmae*eye(o)*randn(o,1)); % Generate white noise sequence with covariance sigma^2*I
    phi_DM(:,k) = H*u(:,k-1);
    eps_true(:,k) = phiSim(:,k)-phi_DM(:,k);
    eps_true(:,k) = eps_true(:,k) - mean(eps_true(:,k));
    s(:,k) = G*(eps_true(:,k)) + e;
    if k == 2
        eps_est(:,2) = B*[zeros(m2,1) ; u(:,2)] + K*s(:,2);
    end
    y(:,k) = K*s(:,k) + (A-K*G)*eps_est(:,k);
    u(:,k) = inv(H)*((A-K*G)*eps_est(:,k) + A*H*u(:,k-1) + K*s(:,k));
    
    eps_est(:,k+1) = (A-K*G)*eps_est(:,k) + B*[u(:,k-1) ; u(:,k)] + K*s(:,k);
    phi_est(:,k) = eps_est(:,k) + H*u(:,k-1);
end
% From the computed matrices the variance is computed. Since the
% variance over N time points is desired, the variance needs to be
% computed per row, or using the var-function with a transposed
% matrix, such that a variance of 1x(m^2) is obtained
VAF = vaf(phiSim,phi_est);
var_eps = var((eps_true-eps_est(:,1:end-1)),0,2);


