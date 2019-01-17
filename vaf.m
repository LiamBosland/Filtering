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

% NUM = (1/N)*norm(y-y_est)^2;
% DEN = (1/N)*norm(y)^2;
% VAF1 = max(0,1-(NUM/DEN))*100;

for k = 1 : N
    yzm(:,k) = y(:,k) - mean(y(:,k));
    y_estzm(:,k) = y_est(:,k) - mean(y_est(:,k));
    dest(:,k) = yzm(:,k) - y_estzm(:,k);
    sub1(1,k) = dest(:,k)'*dest(:,k);
    sub2(1,k) = yzm(:,k)'*yzm(:,k);
end
VAF = max(0,1-(sum(sub1)/sum(sub2)))*100;


