function [var_eps] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim)
% This MATLAB routine closes the loop and computes the variance of the
% residual wavefront var_eps, from system matrices A, G, H, covariance
% matrices Cphi0, sigmae*I, Cw, Kalman-gain matrix K and the wavefront data
% phiSim

% Define some shorthand matrices
B1 = [A*H -H];
F = -B1;
B2 = [(A*H-K*G*H) -H];
%y(k) = K*s(k) + (A-K*G)*eps_est(k);
%% Using the Kalman-filter the loop is closed, computing an estimate for epsilon
% First, it is important to check whether  matrix F is full rank, on which
% the estimation of control action u depends
[U,S,Vt] = svd(F);
Lf = length(nonzeros(round(diag(S),10))); % Find the number of nonzero elements of the singular values of F rounded up to 10 decimal points
if Lf == size(F,1) % F is full row-rank
    FR = 1;
else % F is not full row-rank
    FR = 0;
    U1 = U(:,1:L);
    S1 = S(1:L,1:L);
    V1 = Vt(1:L,:)';
end

% Next, for each set of data from phiSim, the closed-loop is computed
[var_eps_true,var_eps_est] = deal(zeros(length(phiSim),size(phiSim{1,1},1)));
for j = 1 : length(phiSim)
    phi = phiSim{1,j};
    [m2,N] = size(phi);     % The size of phi(k) is (m^2)x1 (real), m = p+1
    [y,u,eps_est,phi_DM] = deal(zeros(m2,N));
    s = zeros(size(G,1),N); % The size of s(k) is 2(p^2)x1 (real)
    for k = 2 : N
        phi_DM(:,k) = H*u(:,k-1);
        s(:,k) = G*(phi(:,k) - phi_DM(:,k));
        y(:,k) = K*s(:,k) + (A-K*G)*eps_est(:,k);
        if FR == 1 % F is full row-rank
            X = pinv(F)*y(:,k);
            u(:,k-1) = X(1:m2);
            u(:,k) = X(m2+1:end);
        elseif FR == 0 % F is not full row-rank
            X = V1*(S1\(U1'))*(U1')*y(:,k);
            u(:,k-1) = X(1:m2);
            u(:,k) = X(m2+1:end);
        end
        if k == 2
            eps_est(:,2) = B2*[zeros(m2,1) ; u(:,2)] + K*G*phi(:,2);
        end
        eps_est(:,k+1) = (A-K*G)*eps_est(:,k) + B2*[u(:,k-1) ; u(:,k)] + K*G*phi(:,k);
    end
    % From the computed matrices the variance is computed. Since the
    % variance over N time points is desired, the variance needs to be
    % computed per row, or using the var-function with a transposed
    % phi-matrix, such that a variance of 1x(m^2) is obtained
    var_eps_true(j,:) = var((phi-phi_DM),0,2);
    var_eps_est(j,:) = var((eps_est(:,1:end-1)),0,2);
end
var_eps = [var_eps_true ; var_eps_est];
end

