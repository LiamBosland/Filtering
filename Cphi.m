function [Cphi_tau] = Cphi(phi,tau)
% Cphi builds the covariance matrix of input phi and time-shift tau,
% according to the approximation:
%
% C_phi(tau) = (1 / (N - tau))*sum_{i=tau+1}^{N}(phi(i)*(phi(i-tau)')
% Input phi is a column vector. If a row vector is supplied, it is
% transposed to a column vector.

% If only one input is given, thus tau is zero
if nargin == 1
    tau = 0;
end
[N,lphi] = size(phi);
% If phi is in the form of vectors, but multiple (matrix). Vector input
% is required.
if N ~= 1 && lphi > 1
    error('The input phi is not a vector. Either vectorize phi or supply approriate input');
end
    
% If a row vector is supplied, transpose it
if lphi > N && N == 1
    phi = phi';
end
N = length(phi);
Cphi_tau = (1/(N-tau))*(phi(tau+1:end)*(phi(1:end-tau)'));
end

