
%% Final Assignment SC42025 Filtering and Identification
% Written by Simon Stouten (4195981) and Liam Bosland(999999)
%
% The contents to this assignment are:
%
%
%
%
% Practical Assignment 2018 - 2019
% Turbulence Modeling for Adaptive Optics
clc;
clearvars;
close all;
tic;
%% Given data
load systemMatrices.mat
load turbulenceData.mat

%% Model 1: Random-walk Model

%% Model 2: Vector Auto-Regressive Model of Order 1
% To formulate the VAR model of order 1 a number of steps have been taken
% which are not discussed in this MATLAB code. For the full explanation and
% elaboration to this model the reader is referred to the accompanying
% report (pXX)

% TEST: Take the first dataset of PhiIdent to test the model
for i = 1 : 1 %length(phiIdent)
    phi = phiIdent{1,i};
    % Since full acces to phi is available, it can be used to approximate
    % the covariance matrices Cphi
    Cphi0 = Cphi(phi);
    Cphi1 = Cphi(phi);
    [A,Cw,K] = computeKalmanAR(Cphi0,Cphi1,G,sigmae);
end

[var_eps] = AOloopAR(G,H,Cphi0,sigmae,A,Cw,K,phiSim);
toc;