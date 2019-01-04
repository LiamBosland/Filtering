function [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae)
% This MATLAB routine computes system state matrix A, the covariance matrix
% Cw of noise signal w(k) and Kalman gain K, from the covariance matrix of
% phi, Cphi, at time instances 0 and 1, output state-matrix C(= G) and the
% variance of the error signal e(k), denoted by sigmae*I

% The first computation is finding A using a least-squares estimate
% Before this can be done however, it is important to check if matrix Chpi0
% is full rank:
[U,S,V] = svd(Cphi0);
L = length(nonzeros(round(diag(S),10))); % Find the number of nonzero elements of the singular values of Cphi0 rounded up to 10 decimal points
if L == size(Cphi0,1) && L == size(Cphi0,2) % Cphi0 is full rank
    At = pinv((Cphi0'))*(Cphi1');
    A = At';
else % Cphi0 is not full rank
    S1 = S(1:L,1:L);
    U1 = U(:,1:L);
    V1 = V(1:L,:)';
    At = V1*inv(S1)*(U1')*Cphi1;
    A = At';
end
% Next is the computation of Cw:
Cw = Cphi0 - (A*Cphi0*(A'));
if issymmetric(Cw) == 0 % Matrix Cw is NOT symmetric, thus we force it to be symmetric
    Cws = 0.5+(Cw + Cw'); % Note that this procedure may introduce errors
    Cw = Cws;             % if original Cw is far from symmetric
elseif issymetric(Cw) == 1 % Matrix Cw is symmetric, which is desired
end

% Before the Kalman-gain matrix can be computed, the two conditions to the
% Kalman filter must be tested:
% -------------------------------------------------------------------------------------------------------------
% 1) The pair (A,G) must be observable.
W = obsv(A,G);
% W must be full-rank in order for the pair to be observable, this can be
% checked using SVD again
[~,Sw,~] = svd(W);
Lw = length(nonzeros(round(diag(Sw),10))); % Find the number of nonzero elements of the singular values of W rounded up to 10 decimal points
if Lw == size(W,2) % W is full column-rank
else % W is not full column-rank
   error('The pair (A,G) is not observable, provide different system matrices');  
end

% 2) The pair (A,sqrt(Q)) must be reachable. If its controllability
% matrix full rank this implies reachability. For this case, Q = Cw
Ctr = ctrb(A,sqrt(Cw));
% Ctr must be full-rank in order for the pair to be controllabel thus 
% reachable, this can be checked using SVD again
[~,Sc,~] = svd(Ctr);
Lc = length(nonzeros(round(diag(Sc),10))); % Find the number of nonzero elements of the singular values of Ctr rounded up to 10 decimal points
if Lc == size(Ctr,1) % Ctr is full row-rank
else % Ctr is not full row-rank
   error('The pair (A,sqrt(Cw)) is not reachable, provide different system matrices');  
end
% -------------------------------------------------------------------------------------------------------------

% For DARE to come to a solution, Cw must be a positive (semi-)definite
% matrix. Positive definiteness can be checked using the CHOL-funtion:
[~,psd] = chol(Cw);
if psd == 0 % Matrix Cw is positive definite
elseif psd > 0 % Matrix Cw is NOT positive definite
    disp('Computed matrix Cw is not positive definite!');
end

% If both conditions hold, the Kalman-gain matrix can be computed using the DARE:
Q = abs(Cw);
R = sigmae*eye(size((G'),2));
[P,~,Kt,~] = dare(A',G',Q,R);
K2 = (A*P*(G'))*inv(R+G*P*(G'));
K = Kt';
end

