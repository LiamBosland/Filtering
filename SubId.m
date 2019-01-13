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
% Construct past and future Hankel matrices of output signal s(k):
Nh = size(s_id,2)/2;
Sp = Hankel(s_id,1,r,Nh);
Sf = Hankel(s_id,r,r,Nh);
[mp,pp] = size(Sp); [mf,pf] = size(Sf);
% Check if the Hankel matrices have equal size:
if mp == mf && pp == pf
else
    error('Hankel matrices did not have equal size, check dimensions');
end
% Compute the RQ-factorization using the qr-function
[R,~] = qr([Sp;Sf]);
R11 = R(1:mp,1:pp);
R21 = R(mp+1:end,1:pp);
%R22 = R(mp+1:end,pp+1:end);
Y = R21*inv(R11)*Sp;
% In order for the RQ-factorization to give a sensible estimation, we
% require R21*inv(R11)*Sp to have rank n which is checked using the SVD:
[Uy,Sy,Vy] = svd(Y);
Ly = length(nonzeros(round(diag(Sy),10))); % Find the number of nonzero elements of the singular values of Y rounded up to 10 decimal points
if Ly == n
elseif Ly ~= n
    error('Rank condition of matrix R21*inv(R11)*Sp does not hold, observability matrix O cannot be computed');
end
% The observability matrix can be found in matrix U (up to a similarity
% transformation):
Cs = Uy(1:l,:);

% Matrix A_s is found through the solution of an least-squares problem:
% min_A e^T*e, s.t. e = y - FA
y = Uy(l+1:end,:);
F = Uy(1:l*(r-1),:);
% Check if matrix F is full-rank, on which the solution of A depends:
[Uf,Sf,Vf] = svd(F);
Lf = length(nonzeros(round(diag(Sf),10))); % Find the number of nonzero elements of the singular values of F rounded up to 10 decimal points
if Lf == size(F,2) % F is full column-rank (rank n)
    As = pinv((F))*(y);
else % F is not full column-rank
    S1 = Sf(1:Lf,1:Lf);
    U1 = Uf(:,1:Lf);
    V1 = Vf(1:Lf,:)';
    As = V1*inv(S1)*(U1')*y;
end

% Finally an estimate for state future Hankel matrix X_(r,N) is computed:
Xf_est = (Sy)^(1/2)*Vy;

end

