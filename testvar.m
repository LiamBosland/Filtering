clc; clearvars; close all;

n = 5000;
N = 10;
sigma = 1;
for k = 1 : N
    e1 = sigma*eye(n)*randn(n,1);
    e2 = sigma*eye(n)*wgn(n,1,1);
    m(k,:) = [mean(e1) mean(e2) var(e1) var(e2)];
end

M = mean(m)