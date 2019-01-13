function [Y] = Hankel(y,k,s,N)
% Hankel(y,k,s,N) builds a Hankel matrix using vectorized input y, the 
% starting index k the amount of block-rows s and the amount of 
% block-columns N

% Default case, y is supplied and a square Hankel matrix is to be built
if nargin == 1
    k = 1;
    s = floor(length(y)/2);
    N = ceil(length(y)/2);
% The case where i is omitted, this the default is i = 1
elseif nargin == 3
    st = k;
    N = s;
    s = st;
    k = 1;
% The case where all information is supplied, this a custom Hankel matrix
% is built
end

% Checking if the inputs are of appropriate dimensions/criteria
[Ny,ly] = size(y);
if ly > 1 && Ny > 1
    error('The input y is not a vector but a matrix. Please vectorize or supply a different input');
elseif ly > 1 && Ny == 1
    error('The input is a row vector, please supply a column vector');
end
% MATLAB cannot take zero indices, thus all k = 0 inputs are set to k = 1
if s == 0
    s = 1;
end
if N+s > Ny
    error('The Hankel dimensions are too large for input y, either provide more data or smaller Hankel parameters s and N');
end

c = 0;
Y = zeros(s,N);
    for i = k : s
        for j = k : N
            Y(i,j) = y(j+c);
        end
        c = c + 1;
    end
end