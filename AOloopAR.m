function [var_eps] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim)
% This MATLAB routine closes the loop and computes the variance of the
% residual wavefront var_eps, from system matrices A, G, H, covariance
% matrices Cphi0, sigmae*I, Cw, Kalman-gain matrix K and the wavefront data
% phiSim

% Define some shorthand matrices
B1 = [A*H -H];
F = -B1;
B2 = [(A*H-K*H) -H]; 
%y(k) = K*s(k) + (A-K*G)*eps_est(k);
%% Using the Kalman-filter the loop is closed, computing an estimate for epsilon
% First, it is important to check whether  matrix F is full rank, on which
% the estimation of control action u depends
[U,S,Vt] = svd(F);
Lb = length(nonzeros(round(diag(Sb),10))); % Find the number of nonzero elements of the singular values of F rounded up to 10 decimal points
if Lb == size(F,1) % F is full row-rank
    FR = 1;
else % F is not full row-rank
    FR = 0;
    U1 = U(:,1:L);
    S1 = S(1:L,1:L);
    V1 = Vt(1:L,:)';  
end
for j = 1 : 1 %length(phiSim)
    phi = phiSim{1,j};
    N = length(phi);
    [s,y,u,eps_est,phi_DM] = deal(zeros(N+1,1));
    for k = 2 : N+1
        if FR == 1 % F is full row-rank
            u(k) = pinv(F)*y(k);
        elseif FR == 0 % F is not full row-rank
            u(k) = V1*(S1\(U1'))*(U1')*y(k);
        end
        if k == 2
            disp('K2');
            eps_est(2) = B2*[0 ; u(2)] + K*phi(2);
        end
        phi_DM(k) = H*u(k-1);
        s(k) = G*(phi(k) - phi_DM(k));
        y(k) = K*s(k) + (A-K*G)*eps_est(k);
        eps_est(k+1) = (A-K*G)*eps_est(k) + B*[u(k-1) ; u(k)] + K*s(k);
    end
    var_eps = [var(eps_est) var(phi-phi_DM)];
end

end

