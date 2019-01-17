function [VAF] = vaf(y,y_est)
% This function returns the Variance Accounted For (VAF) of estimate y_est,
% by comparing it to measured signal y

[ny,Ny] = size(y); [n,N] = size(y_est);
if ny == n && Ny == N
else
    error('The provided inputs do not have the same size.');
end
if n > N % y_est has more rows than columns, so it will be transposed
    y = y';
end

NUM = (1/N)*norm(y-y_est)^2;
DEN = (1/N)*norm(y)^2;
VAF = max(0,1-(NUM/DEN));


