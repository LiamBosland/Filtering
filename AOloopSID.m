function [var_eps,VAF] = AOloopSID(G,H,As,Cs,Ks,sigmae,phiSim)
% This MATLAB routine closes the loop for the Subspace Identification model
% and computes the variance of residual wavefront, using system matrices
% As, Cs, Ks, G, H, covariance matric sigmae^2*I and simulated data phiSim

% Define some shorthards
o = size(G,1);  % The amount of outputs per time-step k (equal to 2*p^2)
Cphi0 = Cphi(phiSim);
Gamma = Cphi0*G'*(G*Cphi0*G' + sigmae^2*eye(o));
G1 = Gamma*Cs;
G2 = Ks*pinv(Gamma);
A_star = As*pinv(G1) - G2;
B = [A_star*H -pinv(G1)*H];

% Next, the closed-loop is computed for the provided phiSim
[m2,N] = size(phiSim);     % The size of phi(k) is (m^2)x1 (real), m = p+1
[u,eps_est,eps_true,phi_DM,phi_est] = deal(zeros(m2,N));
sc = zeros(o,N);
so = [(G*phiSim(:,1) + (sigmae^2*eye(o)*randn(o,1))) zeros(o,N-1)];
for k = 2 : N
    e = sigmae^2*eye(o)*randn(o,1);
    so(:,k) = G*phiSim(:,k) + e;
    if k == 2
        eps_est(:,2) = G1*(B*[zeros(m2,1) ; u(:,2)]+ Ks*so(:,2));
    end
    phi_DM(:,k) = H*u(:,k-1);
    eps_true(:,k) = phiSim(:,k) - phi_DM(:,k);
    eps_true(:,k) = eps_true(:,k) - mean(eps_true(:,k));
    sc(:,k) = pinv(Gamma)*(eps_est(:,k) + H*u(:,k));
    u(:,k) = inv(H)*G1*(A_star*(eps_est(:,k) + H*u(:,k-1)) + Ks*so(:,k));
    eps_est(:,k+1) = G1*(A_star*eps_est(:,k) + ...
        B*[u(:,k-1) ; u(:,k)] + Ks*so(:,k));
    phi_est(:,k) = eps_est(:,k) + H*u(:,k-1);
end
VAF = vaf(phi,phi_est)
var_eps = var((eps_true-eps_est(:,1:end-1)),0,2);
