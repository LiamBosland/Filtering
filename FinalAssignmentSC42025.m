%% Final Assignment SC42025 Filtering and Identification
% Written by Simon Stouten (4195981) and Liam Bosland(4216377)
%
% The contents to this assignment are:
% AOloop_nocontrol.m
% AOloopAR.m
% AOloopRW.m
% AOloopRWSlopes.m
% AOloopSID.m
% computeKalmanAR.m
% Cphi.m
% Hankel.m
% matstable.m
% Rank.m
% SubId.m
% vaf.m
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
ns = length(phiSim);    % The number of samples available for analysis
%% Model 1: Random-walk Model
disp('Evaluating the Random-walk Model');
tic;
VAF_RW = zeros(ns,1);
for i = 1 : ns
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1);
    [var_eps_RW(i,:),VAF_RW(i,1)] = AOloopRW(G,H,Cphi0,sigmae,phisim);
end
%% Model 2: Vector Auto-Regressive Model of Order 1
disp('Evaluating the VAR1 Model');
% Define zero-matrices
[stable_AR,VAF_AR] = deal(zeros(ns,1));
for i = 1 : ns
    phisim = phiSim{1,i};
    % Since full acces to phi is available, it can be used to approximate
    % the covariance matrices Cphi
    Cphi0 = Cphi(phisim);
    Cphi0 = 0.5*(Cphi0 + Cphi0');
    Cphi1 = Cphi(phisim,1);
    [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae);
    stable_AR(i,1) = matstable(A-K*G);
    [var_eps_AR(i,:),VAF_AR(i,1)] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phisim);
end
%% Determining order range of the system
fig1 = figure;
for k = 1 : ns
    phi = phiIdent{1,k};
    Y = Hankel(phi,1,10,2500);
    [~,S,~] = svd(Y);
    semilogy(S,'b*'); grid on; hold on; legend on;
    xlabel('Singular value no.'); ylabel('Intensity singular value');
    title('Singular values of phiIdent');
end
saveas(fig1,'SingularValuesPhi','png');
clear fig1
%% Model 3: Subspace Identification
disp('Evaluating the Subspace Identification routine');
% For identification 2/3 of the data is used, the remaining 1/3 is used for
% the validation process, using cross-validation.
n_id = 4;
n = linspace(100,250,n_id);
VAR = cell(length(phiIdent),length(n));
stable_SI = zeros(length(phiIdent),length(n));
for j = 1 : n_id
    for i = 1 : ns
        fprintf('SubID: order = %0i, Sample %0i',n(j),i); fprintf('\n');
        phi = phiIdent{1,i};
        sysorder = n(j); N = size(phi,2); o = size(G,1); % Define some dimensions
        r = 10;  % Since a MIMO system is to be identified, we require 2p^2*r > n
        N_id = 3500; N_val = N-N_id; % Approximately 2/3, 1/3
        % Simulate open-loop measurements so(k):
        s_id = zeros(o,N);
        for k = 1 : N
            e = (sigmae^2*eye(o)*randn(o,1)); % Generate white noise sequence with covariance sigma^2*I
            s_id(:,k) = G*phi(:,k) + e;
        end
        [As,Cs,Ks] = SubId(s_id,N_id,N_val,r,sysorder);
        stable_SI(i,j) = matstable(As-Ks*Cs);
        [var_eps(i,j),VAF_SI(i,j)] = AOloopSID(G,H,As,Cs,Ks,sigmae,phi);
    end
end

%% Model 4: Random-walk Model Residual Slopes
disp('Evaluating the Random-walk with Residual Slopes model');
VAF_SL = ones(ns,1);
for i = 1 : ns
    phisim = phiSim{1,i};
    Cphi0 = Cphi(phisim);
    Cphi1 = Cphi(phisim,1);
    [var_eps_RSL(i,:),VAF_SL(i,1)] = AOloopRWSlopes(G,H,Cphi0,sigmae,phisim);
end
toc;
%% Model 5: No control
disp('Evaluating the no control model');
VAF_NC = ones(ns,1);
for i = 1 : ns
    phisim = phiSim{1,i};
    [sigma(i,1),var_NC(i,:)] = AOloop_nocontrol(phisim,sigmae,H,G);
end
toc;

%% Evaluating Model Performances
for j = 1 : n_id
    maxsi(j) = max(VAF_SI(:,j));
    meansi(j) = mean(VAF_SI(:,j));
    [~,I] = min(mean(cell2mat(VAR)));
end
var_SI_best = cell2mat(VAR(:,I));
[~,I] = max(meansi);
VAF_SI_best(:,1) = VAF_SI(:,I);

% Finding the best variance among the 20 samples, by evaluating maximum
% mean
% meanvar_NC = mean(var_NC,1); [~,I] = min(meanvar_NC); 
% varNC = var_NC(:,I);
% meanvar_RW = mean(var_eps_RW,1); [~,I] = min(meanvar_RW); 
% var_RW = var_eps_RW(:,I);
% meanvar_AR = mean(var_eps_AR,1); [~,I] = min(meanvar_AR); 
% var_AR = var_eps_AR(:,I);
% meanvar_SI = mean(var_eps_SI,1); [~,I] = min(meanvar_SI); 
% var_SI = var_eps_SI(:,I);
% meanvar_RSL = mean(var_eps_RSL,1); [~,I] = min(meanvar_RSL); 
% var_RSL = var_eps_RSL(:,I);

fig2 = figure;
plot(n,meansi,'r*'); grid on;
xlabel('Sample no.'); ylabel('Variance Accounted For [%]');
title('VAF Subspace Identification vs. system order n');
saveas(fig2,'MeansVAFSI','png');
clear fig2

fig3 = figure;
hold on; grid on;
bar(1:20,[VAF_RW VAF_AR VAF_SI_best VAF_SL]);
legend('Random-walk','VAR1','SID','Residual Slopes');
xlabel('Sample no.'); ylabel('Variance Accounted For [%]');
title('VAF of given data samples');
saveas(fig3,'VAFComparison','png');
clear fig3

fig4 = figure;
hold on; grid on
plot((1:20),sigma,'k*',(1:20),var_eps_RW,'ro',(1:20),var_eps_AR,'gx',...
    (1:20),var_SI_best,'md',(1:20),var_eps_RSL,'r*');
legend('No control','Random-walk','VAR1','SID','Residual Slopes');
xlabel('Sample no.'); ylabel('Variance estimated epsilon');
title('Variance estimated residual wavefront');
saveas(fig4,'VARComparison','png');
clear fig4
%% Identifying unobservable mode
% It is given that two unobservable modes exist in the system, which can
% be identified through matrix G. One mode is given to be phi(k) =
% ones(49,1), related to the piston not being measured. The other mode will
% likely be an unobservable wavefront. First, the null-space of G is
% computed
NullG = null(G);

% Since the first move must be related to the piston, we must search for a
% transformation vector T such that NullG(:,1).*T = ones(49,1). Since we
% can identify two distinct values in NullG(:,1), T will be a repetition of
% 1 devided by these factors:
a1 = 1/NullG(1,1); a2 = 1/NullG(2,1);
T = [repmat([a1;a2],[24,1]) ; a1];
% Check if NullG(:,1).*T = ones(49,1):
unob1 = NullG(:,1).*T;
% Now the second column of NullG can be mapped using T and thus the
% unobservable mode can be reconstructed:
unob2 = NullG(:,2).*T;
fig5 = figure;
imagesc(reshape(unob2,[7,7]));
saveas(fig5,'UnobservableMode21','png');
clear fig5
N = null([G;ones(1,49)]);
fig6 = figure;
imagesc(reshape(N,[7,7]));
saveas(fig6,'UnobservableMode22','png');
clear fig6
%% Saving Workspace Data
save('WorkspaceFiltering.mat');
pause;