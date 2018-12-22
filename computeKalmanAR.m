function [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae)
% This MATLAB routine computes system state matrix A, the covariance matrix
% Cw of noise signal w(k) and Kalman gain K, from the covariance matrix of
% phi, Cphi, at time instances 0 and 1, output state-matrix C(= G) and the
% variance of the error signal e(k), denoted by sigmae

% The first computation is finding A using a least-squares estimate:
At = pinv((Cphi0'))*(Cphi1');
A = At';

% Next is the computation of Cw:
Cw = Cphi0 - (A*Cphi0*(A'));

% Finally the Kalman gain-matrix can be computed using the DARE:
[~,~,Kt] = dare(A',G',Cw,sigmae*eye(size(G,2)));
K = Kt';
end

