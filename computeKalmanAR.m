function [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae)
% This MATLAB routine computes system state matrix A, the covariance matrix
% Cw of noise signal w(k) and Kalman gain K, from the covariance matrix of
% phi, Cphi, at time instances 0 and 1, output state-matrix C(= G) and the
% variance of the error signal e(k), denoted by sigmae*I

% The first computation is finding A using a least-squares estimate
% Before this can be done however, it is important to check if matrix Chpi0
% is full rank, which is performend by Rank.m:
[rnk,FR,U,S,Vt,~] = Rank(Cphi0,0,true);
if FR == 1 % Cphi0 is full rank
    At = pinv((Cphi0'))*(Cphi1');
    A = At';
elseif FR == 0
    S1 = S(1:rnk,1:rnk);
    U1 = U(:,1:rnk);
    V1 = Vt(1:rnk,:)';
    At = V1*inv(S1)*(U1')*Cphi1;
    A = At';
end
% Next is the computation of Cw:
Cw = Cphi0 - (A*Cphi0*(A'));
if issymmetric(Cw) == 0 % Matrix Cw is NOT symmetric, thus we force it to be symmetric
    Cws = 0.5*(Cw + Cw'); % Note that this procedure may introduce errors
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
[~,FR,~,~,~,~] = Rank(W);
if FR == 1 % W is full column-rank
elseif FR == 0 % W is not full column-rank
   error('The pair (A,G) is not observable, provide different system matrices');  
end

% 2) The pair (A,sqrt(Q)) must be reachable. If its controllability
% matrix full rank this implies reachability. For this case, Q = Cw
Ctr = ctrb(A,sqrt(Cw));
% Ctr must be full-rank in order for the pair to be controllabel thus 
% reachable, this can be checked using SVD again
[~,FR,~,~,~,~] = Rank(Ctr);
if FR == 1 % Ctr is full row-rank
elseif FR == 0 % Ctr is not full row-rank
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
Q = Cw;
R = sigmae*eye(size((G'),2));
[~,~,Kt,~] = dare(A',G',Q,R);
K = Kt';
end

