function [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae)
% This MATLAB routine computes system state matrix A, the covariance matrix
% Cw of noise signal w(k) and Kalman gain K, from the covariance matrix of
% phi, Cphi, at time instances 0 and 1, output state-matrix C(= G) and the
% variance of the error signal e(k), denoted by sigmae

% The first computation is finding A using a least-squares estimate
% Before this can be done however, it is important to see if matrix Chpi0
% is full rank:
[U,S,V] = svd(Cphi0);
L = length(nonzeros(round(diag(S),10))); % Find the number of nonzero elements of the singular values of Cphi0 rounded up to 10 decimal points
if L == size(Cphi0,1) && L == size(Cphi0,2) % Cphi0 is full rank
    At = pinv((Cphi0'))*(Cphi1');
    A = At';
else % Cphi0 is not full rank
    S = S(1:L,1:L);
    U1 = U(:,1:L);
    V1 = V(1:L,:)';
    A = V1*(S^-1)*(U1')*Cphi1;
end
% Next is the computation of Cw:
Cw = Cphi0 - (A*Cphi0*(A'));
issymmetric(Cw)
% Finally the Kalman gain-matrix can be computed using the DARE:
[~,~,Kt] = dare(A',G',Cw,sigmae*eye(size(G,2)));
K = Kt';
end

