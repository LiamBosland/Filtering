function [var_eps] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim)
% This MATLAB routine closes the loop and computes the variance of the
% residual wavefront var_eps, from system matrices A, G, H, covariance
% matrices Cphi0, sigmae*I, Cw, Kalman-gain matrix K and the wavefront data
% phiSim

% Define some shorthands
B1 = [A*H -H];
F = -B1;
B2 = [(A*H-K*G*H) -H];
o = size(G,1);  % The amount of outputs per time-step k (equal to 2*p^2)
%y(k) = K*s(k) + (A-K*G)*eps_est(k);
%% Using the Kalman-filter the loop is closed, computing an estimate for epsilon
% First, it is important to check whether  matrix F is full rank, on which
% the estimation of control action u depends
% [U,S,Vt] = svd(F);
% Lf = length(nonzeros(round(diag(S),10))); % Find the number of nonzero elements of the singular values of F rounded up to 10 decimal points
% if Lf == size(F,1) % F is  full row-rank
%     FR = 1;
% else % F is not full row-rank
%     FR = 0;
%     U1 = U(:,1:L);
%     S1 = S(1:L,1:L);
%     V1 = Vt(1:L,:)';
% end

% Next, the closed-loop is computed for the provided phiSim
[m2,N] = size(phiSim);     % The size of phi(k) is (m^2)x1 (real), m = p+1
[y,u,eps_est,eps_true,phi_DM] = deal(zeros(m2,N));
s = [(G*phiSim(:,1) + (sigmae*eye(o)*randn(o,1))) zeros(o,N-1)]; % The size of s(k) is 2(p^2)x1 (real)
for k = 2 : N
    e = (sigmae*eye(o)*randn(o,1)); % Generate white noise sequence with covariance sigma*I
    phi_DM(:,k) = H*u(:,k-1);
    eps_true(:,k) = phiSim(:,k)-phi_DM(:,k);
    eps_true(:,k) = eps_true(:,k) - mean(eps_true(:,k));
    s(:,k) = G*(eps_true(:,k)) + e;
    y(:,k) = K*s(:,k) + (A-K*G)*eps_est(:,k);
        
    %----------------- u1 ----------------------------
%     if FR == 1 % F is full row-rank
%         X = pinv(F)*y(:,k);
%         u1(:,k-1) = X(1:m2);
%         u1(:,k) = X(m2+1:end);
%     elseif FR == 0 % F is not full row-rank
%         X = V1*inv(S1)*(U1')*y(:,k);
%         u1(:,k-1) = X(1:m2);
%         u1(:,k) = X(m2+1:end);
%     end
    if k == 2
        eps_est(:,2) = B2*[zeros(m2,1) ; u(:,2)] + K*s(:,2);
    end
    %----------------- u2 ----------------------------
    u(:,k) = inv(H)*((A-K*G)*eps_est(:,k) + A*H*u(:,k-1) + K*s(:,k));    
    eps_est(:,k+1) = (A-K*G)*eps_est(:,k) + B2*[u(:,k-1) ; u(:,k)] + K*s(:,k);
end
% % From the computed matrices the variance is computed. Since the
% % variance over N time points is desired, the variance needs to be
% % computed per row, or using the var-function with a transposed
% % matrix, such that a variance of 1x(m^2) is obtained
% var_eps = [var(eps_true,0,2) ; var(eps_est(:,2:end),0,2)];

% From the computed residual wavefronts the variance ||eps_true - eps_est||^2_2
% is computed. The variance is computed using the cov-function, to give the
% variance of a vector/matrix.
var_eps = cov(eps_true'-eps_est(:,1:end-1)');


