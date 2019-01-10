function [Cphi_tau] = Cphi(phi,tau)
% Cphi builds the covariance matrix of input phi and time-shift tau
% according to the approximation:
%
% C_phi(tau) = (1 / (N - tau))*sum_{i=tau+1}^{N}(phi(i)*(phi(i-tau)')
% The input phi is a n x N matrix, containing N time samples of n sensors
% and is assumed to be of non-zero mean. Thus the mean is first subtracted
% from the input phi, columnwise.

% If only one input is given, thus tau is zero
if nargin == 1
    tau = 0;
end
[nphi,N] = size(phi);
% The amount of time samples is smaller than the amount of sensors, this
% either means a bad dataset is provided or that the input is transposed.
% It is transposed again to try and solve the issue.
if nphi > N
    disp('The input has more rows than columns, indicating either a bad dataset or a transposed one.');
    disp('The input is transposed. If the matrix was provided in correct dimensions, check the validity of the time samples');
    phi = phi';
    [nphi,N] = size(phi);
end

% Removing the mean from input phi, columnwise.
mu = mean(phi,1);
phi_zm = phi - mu;
p1 = phi(:,tau+1:end);
p2 = phi(:,1:N-tau);

% Computing the Covariance matrix Cphi
Cphi_tau = (1/(N-tau))*(p1*p2');

