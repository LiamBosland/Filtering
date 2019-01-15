function [As,Cs,Ks] = SubId(s_id,N_id,N_val,r,n)
% This MATLAB routine identifies matrices A_s, C_s and K_s from open-loop
% wavefront sensor data s_id, N_id and N_val, the number of samples for 
% identification and validation, respectively. Parameter r is the upper
% bound on system order n. 
% 
% Note that in the original assignment formulation
% and in the book, parameter s is the upper bound to n. For the sake of
% clarity the parameter r has been used for this report instead, since the
% output is stored in the slopes s(k). Thus to avoid confusion between the
% upper bound s and slopes s(k), r was chosen as representation for s

% Check if the supplied dimensions of s_id are correct, it should be 
% [2p^2 x N]
if size(s_id,1) > size(s_id,2) % s_id has more entries in depth than width
    s_id = s_id';
end
l = size(s_id,1);   % stores the amount of outputs at time instant k

% First, the measurements s_id are split up into a identification sequence
% si(k) and a validation sequence s_val(k)
si = s_id(:,1:N_id); 
s_val = s_id(:,N_id+1:end);

% Construct past and future Hankel matrices of output signal s(k):
Nh = size(si,2)/2;
if isinteger(Nh) == 0
    Nh = floor(Nh);
end
Sp = Hankel(si,1,r,Nh);
Sf = Hankel(si,r,r,Nh);
[mp,pp] = size(Sp); [mf,pf] = size(Sf);
% Check if the Hankel matrices have equal size:
if mp == mf && pp == pf
else
    error('Hankel matrices did not have equal size, check dimensions');
end
p = pp/2;
if mod(p,1) ~= 0 % p has decimals
    error('Dimensions do not line up nicely, please check N_id and N_val again');
end
% Compute the RQ-factorization using the qr-function
R = triu(qr([Sp;Sf]'))';
% A = R(1:mp,1:mp); B = R(1:mp,mp+1:end);
% C = R(mp+1:end,1:mp); D = R(mp+1:end,mp+1:end);
% [mean(mean(A,2)); mean(mean(B,2)); mean(mean(C,2)); mean(mean(D,2))]

R11 = R(1:mp,1:p);
R21 = R(mp+1:end,1:p);
%R22 = R(mp+1:end,pp+1:end);
Y = R21*pinv(R11)*Sp;
% In order for the RQ-factorization to give a sensible estimation, we
% require R21*inv(R11)*Sp to have rank n which is checked using the SVD:
[Uy,Sy,Vy] = svd(Y);
% Ly = length(nonzeros(round(diag(Sy),10))); % Find the number of nonzero elements of the singular values of Y rounded up to 10 decimal points
% if Ly == n
% elseif Ly ~= n
%     error('Rank condition of matrix R21*inv(R11)*Sp does not hold, observability matrix O cannot be computed');
% end
% The observability matrix can be found in matrix U (up to a similarity
% transformation):
Cs = Uy(1:l,:);

% Matrix A_s is found through the solution of an least-squares problem:
% min_A e^T*e, s.t. e = y - FA
y = Uy(l+1:end,:);
F = Uy(1:l*(r-1),:);
% Check if matrix F is full-rank, on which the solution of A depends:
[UF,SF,VF] = svd(F);
Lf = length(nonzeros(round(diag(Sf),10))); % Find the number of nonzero elements of the singular values of F rounded up to 10 decimal points
if Lf == size(F,2) % F is full column-rank (rank n)
    As = pinv((F))*(y);
else % F is not full column-rank
    S1 = SF(1:Lf,1:Lf);
    U1 = UF(:,1:Lf);
    V1 = VF(1:Lf,:)';
    As = V1*inv(S1)*(U1')*y;
end

% Finally an estimate for future state sequence X_(r,N) is computed:
Xf_est = (Sy).^(1/2)*Vy;
X_est1 = Xf_est(:,1:end-1); X_est2 = Xf_est(:,2:end); len = size(X_est1,2);

% Next, the model residuals can be computed, giving the Q,S,R covariance
% matrices:
S = Hankel(si,r,1,N_id-r-1);
WV = [X_est2;S] - [As;Cs]*X_est1;
W = WV(1:n,1:len); V = WV(n+1:end,1:len);
[mw,pw] = size(W); [mv,pv] = size(V);
% Check if the matrices have equal size:
if mw == mv && pw == pv
else
    error('Matrices did not have equal size, check dimensions');
end
C = cov(W,V);
Q_est = C(1:mw,1:mw); S_est = C(1:mw,mw+1:end); 
St_est = C(mw+1:end,1:mw); R_est = C(mw+1:end,mw+1:end);
% These matrices should form a symmetric block-matrix [Q S; S^T R], which
% can be checked:
if S_est == transpose(St_est)
else
    error('Computed covariances from W and V are not symmetric');
end
% Matrices Q,S,R can be used to compute Kalman-gain matrix K, using the
% DARE-equation:
[~,~,Kt] = dare(As',Cs',Q_est,R_est);
Ks = Kt';
end

