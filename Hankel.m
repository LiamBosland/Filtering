function [Y] = Hankel(y,k,s,N)
% Hankel(y,k,s,N) builds a Hankel matrix using vectorized input y, the
% starting index k the amount of block-rows (of dimension ly) s and the
% amount of block-columns N.

% Default case, y is supplied and a square Hankel matrix is to be built
if nargin == 1
    k = 1;
    s = floor(length(y)/2);
    N = ceil(length(y)/2);
    % The case where k is omitted, the default is k = 1
elseif nargin == 3
    st = k;
    N = s;
    s = st;
    k = 1;
    % The case where all information is supplied, this a custom Hankel matrix
    % is built
end

% Checking if the inputs are of appropriate dimensions/criteria.
% Time-samples are required to be in the width of supplied y.
[ly,Ny] = size(y);
if ly == 1 || Ny == 1 % Supplied y is a vector
    M = s-1;
    if ly > Ny
        y = y'; % Transpose y such that a column vector is obtained
        [ly,Ny] = size(y);
    end
else % Supplied y is a matrix
    M = 0;
    if ly > Ny % Matrix y is deeper than wider (time-samples required in width)
        y = y'; % Transpose y such that time-samples (assumed to be larger than measurement points y(k)) are in the width of y
        [ly,Ny] = size(y);
    end
end
% MATLAB cannot take zero indices, thus all k = 0 inputs are set to k = 1
if s == 0
    s = 1;
end
if N+s > Ny
    error('The Hankel dimensions are too large for input y, either provide more data or smaller Hankel parameters s and N');
end

c = 0;
Y = zeros(s*ly,N);
for i = 1 : s
    a1 = (i-1)*(ly-1) + i; a2 = a1 + (ly-1); % Two counting parameters
    for j = 1 : N
        Y(a1:a2,j) = y(:,k+j+c-1);
    end
    c = c + 1;
end
end