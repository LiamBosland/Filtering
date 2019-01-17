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
tic;
clc
load systemMatrices.mat
load turbulenceData.mat

ns = length(phiSim);
VAF = zeros(1,ns);
for i = 1:ns
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1); 
    
    [eps_est_1,delta_u,s,phi_estRW,VAFRW(i,1)] = AOloopRW(G,H,Cphi0,sigmae,phisim);
end
toc
%% Model 2: Vector Auto-Regressive Model of Order 1
% Define zero-matrices
ns = length(phiSim); % The number of samples available for analysis
stable = zeros(ns,1);
VAF_eps = zeros(size(phiSim{1,1},1),ns);
for i = 1 : ns
    phisim = phiSim{1,i};
    % Since full acces to phi is available, it can be used to approximate
    % the covariance matrices Cphi
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1);
    [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae);
    stable(i) = matstable(A-K*G);
    [var_eps,VAF(i)] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phisim);
end
% Find the highest average VAF among the 20 samples and plot the VAF-vector
% corresponding to this highest average
m = mean(VAF,1);
[VAF_avgmax,I] = max(m);
VAF_best = VAF(:,I);
figure;
grid on;
histogram(VAF_best); xlabel('Variance Accounted For [%]'); ylabel('No. of sensors [-]');
title('VAF of residual wavefront'); legend(strcat('Average: ',num2str(round(VAF_avgmax,1)),' %'));
%% Model 3: Subspace Identification
ns = length(phiIdent);
% For identification 2/3 of the data is used, the remaining 1/3 is used for
% the validation process, using cross-validation.

for i = ns
    phi = phiIdent{1,i};
    n = 50; N = size(phi,2); o = size(G,1); % Define some dimensions
    r = 15;  % For Subspace Identification we require r > n
    N_id = 3500; N_val = N-N_id; % Approximately 2/3, 1/3
    % Simulate open-loop measurements so(k):
    s_id = zeros(o,N);
    for k = 1 : N
        e = (sigmae^2*eye(o)*randn(o,1)); % Generate white noise sequence with covariance sigma^2*I
        s_id(:,k) = G*phi(:,k) + e;
    end
    [As,Cs,Ks] = SubId(s_id,N_id,N_val,r,n);
    [var_eps,VAF(i,:)] = AOloopSID(G,H,As,Cs,Ks,sigmae,phi);
end

%% Model 4: Random-walk Model Residual Slopes
tic;
clc
load systemMatrices.mat
load turbulenceData.mat

ns = length(phiSim);
VAFSL = ones(20,1);
for i = 1:1
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1); 
    
    [eps_est_1,delta_u,s,phi_estSL,VAFSL(i,1)] = AOloopRWSlopes(G,H,Cphi0,sigmae,phisim);
end
toc
%% Make a moving picture for motivation
close all

for Frame = 1:5000
    for k = 1:7
        for i = 1:7
            pic(i,k) = phi_estSL(i+((k-1)*7),Frame);
        end
    end
% subplot(5,9,Frame)
% imagesc(pic);
figure(1)
pause(0.02)
imagesc(pic)
axis off
end

%% Identify nullspace
Null = null([G; ones(49)]');
Zero = G*ones(49);

