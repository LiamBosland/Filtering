%% Final Assignment SC42025 Filtering and Identification
% Written by Simon Stouten (4195981) and Liam Bosland(4216377)
%
% The contents to this assignment are:
%
%
%
%
% Practical Assignment 2018 - 2019
% Turbulence Modeling for Adaptive Optics
%
% To formulate the discussed models a number of steps have been taken
% which may not be fully discussed in this MATLAB code. For the full
% elaborations to this model the reader is referred to the accompanying
% report (pXX)
clc;
clearvars;
close all;
tic;
%% Given data
load systemMatrices.mat
load turbulenceData.mat

%% Model 1: Random-walk Model
clc
clearvars;
load systemMatrices.mat
load turbulenceData.mat

ns = length(phiSim);
e = cell(1,ns);
for i = 1:ns
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1);
    
    [eps(:,:),var_eps,avg_var_eps] = AOloopRW(G,H,Cphi0,sigmae,phisim);
    e{i} = eps;
    var(:,i) = var_eps;
    avg(:,i) = avg_var_eps;
end

%% Model 2: Vector Auto-Regressive Model of Order 1
% Define zero-matrices
ns = length(phiSim); % The number of samples available for analysis
stable = zeros(ns,1);
for i = 1 : ns
    phisim = phiSim{1,i};
    % Since full acces to phi is available, it can be used to approximate
    % the covariance matrices Cphi
    Cphi0 = Cphi(phisim);
    Cphi0 = 0.5*(Cphi0 + Cphi0');
    Cphi1 = Cphi(phisim,1);
    [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae);
    stable_AR(i,1) = matstable(A-K*G);
    [var_eps_AR,VAF_AR(i,1)] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phisim);
end

%% Model 3: Subspace Identification
ns = length(phiIdent);
% For identification 2/3 of the data is used, the remaining 1/3 is used for
% the validation process, using cross-validation.
n = linspace(50,500,10);
VAR = cell(length(phiIdent),length(n));
stable_SI = zeros(length(phiIdent),length(n));
for j = 1 : length(n)
    for i = 1 : ns
        phi = phiIdent{1,i};
        sysorder = n(j); N = size(phi,2); o = size(G,1); % Define some dimensions
        r = 15;  % For Subspace Identification we require r > n
        N_id = 3500; N_val = N-N_id; % Approximately 2/3, 1/3
        % Simulate open-loop measurements so(k):
        s_id = zeros(o,N);
        for k = 1 : N
            e = (sigmae^2*eye(o)*randn(o,1)); % Generate white noise sequence with covariance sigma^2*I
            s_id(:,k) = G*phi(:,k) + e;
        end
        [As,Cs,Ks] = SubId(s_id,N_id,N_val,r,sysorder);
        stable_SI(i,j) = matstable(As-Ks*Cs);
        [var_eps,VAF_SI(i,1)] = AOloopSID(G,H,As,Cs,Ks,sigmae,phi);
        VAR{i,j} = var_eps;
    end
end
disp('end');
%% Model 4: Random-walk Model Residual Slopes
clc
clearvars;
load systemMatrices.mat
load turbulenceData.mat

ns = length(phiSim);
e = cell(1,ns);
res_slopes = cell(1,ns);
for i = 1:ns
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1);
    
    [eps(:,:),var_eps,res_slopes] = AOloopRWSlopes(G,H,Cphi0,sigmae,phisim);
    e{i} = eps;
    var(:,i) = var_eps;
end

rtrn = 0

%% Evaluation of the model performances
% Now that the VAF of each model and each dataset is computed, they can be
% compared.

toc;