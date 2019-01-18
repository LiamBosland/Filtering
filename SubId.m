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
% Compute the RQ-factorization using the qr-function
[~,R] = qr([Sp;Sf]');
R = R';

R11 = R(1:mp,1:mp);
R21 = R(mp+1:2*mp,1:mp);
Y = R21*inv(R11)*Sp;
% From the RQ-factorization state sequence X_(r,N) can be estimated, by
% computing the SVD of Y. Since the n-largest singular values are needed,
% the command svds is used:
[~,Sn,Vn] = svds(Y,n);
% Finally an estimate for future state sequence X_(r,N) is computed:
Xf_est = sqrtm(Sn)*Vn';
X_est1 = Xf_est(:,1:end-1); X_est2 = Xf_est(:,2:end); len = size(X_est1,2);
S = Hankel(si,r,1,Nh-1);

% Construct LSE formulation (y - F*M)
y = [X_est2' S'];
F = X_est1';
% Check if matrix F is full rank or not
[rnk,FR,UF,SF,VF,~] = Rank(F,0,true);
if FR == 1
    M_est = pinv(F)*y;
    M_est = M_est';
elseif FR == 0
    U1 = UF(1:rnk,:);
    S1 = SF(1:rnk,1:rnk);
    V1 = VF(:,1:rnk);
    M_est = V1*inv(S1)*(U1')*y;
    M_est = M_est';
end
% From the estimation of the matrices, M_est, matrices As and Cs can be
% extracted
As = M_est(1:n,1:n);
Cs = M_est(n+1:n+l,1:n);

% Next, the model residuals can be computed, giving the Q,S,R covariance
% matrices:
WV = [X_est2;S] - [As;Cs]*X_est1;
W = WV(1:n,:); V = WV(n+1:end,:); 
[mw,~] = size(W);

C = (1/N_id)*[W;V]*[W' V'];
Q_est = C(1:mw,1:mw); S_est = C(1:mw,mw+1:end); 
St_est = C(mw+1:end,1:mw); R_est = C(mw+1:end,mw+1:end);
% Matrices Q,S,R can be used to compute Kalman-gain matrix K, using the
% DARE-equation:
[~,~,Kt] = dare(As',Cs',Q_est,R_est);
Ks = Kt';
end

