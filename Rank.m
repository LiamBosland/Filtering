function [rnk,FR,U,S,V,Strue] = Rank(A,red,mout,check)
% This is a MATLAB function to compute the rank of matrix A using the
% Singular Value Decomposition. Input argument red can be provided to 
% compute a reduced svd using svds (default is 0). Input red can also be 
% provided a string or char 'n', this is compute the SVD reduced by the 
% computed rank. If input argument mout is provided as true, Rank also 
% returns the matrices U,S,V' of svd(A). If mout is not provided it is set 
% to default: false. Output argument rnk gives the rank of matrix A, while 
% boolean FR returns (if called to) if matrix A is full-rank
% (row or column) or not
%
% To compensate for rounding errors (10^-14 ~ 10^-16) introduced by MATLAB,
% the singular values are first rounded off to 13 decimal places. This is a
% heuristic approach however and will not give a correct result in 100% of
% the uses. Thus caution is to be taken when using this function. To check
% for such errors output Strue can be returned by providing boolean 'check'
% to be true. The default value for 'check' is false.

if nargin < 2 % no additional input arguments are given, thus set to default
    red = 0;
    [mout,check] = deal(false);
end
if nargin <3 % red is provided, but mout and check aren't, thus set to default
    [mout,check] = deal(false);
end
if nargin <4 % Only check is left out, thus set to default
    [check] = deal(false);
end
[m,n] = size(A);
p = min(m,n);   % rank-check parameter, will equal the lowest value of m,n
[U,S,V] = svd(A);
Strue = S;
v = round(diag(S),13);  % vector v is constructed with rounded-off singular values
L = length(nonzeros(v)); % length of nonzero elements in v is computed
if L == p % A is full rank
    FR = true;
elseif L < p % A is not full rank
    FR = false;
elseif L > p
    error('Something went wrong in computing the rank and singular values, please check matrix A');
end
rnk = L;

if ischar(red) == 1 || isstring(n) == 1 % If input red is a character or string
    red = L;
else
end
if red == 0
elseif red ~= 0
    [U,S,V] = svds(A,red);
end
if mout == false && check == false
    [U,S,V,Strue] = deal([]);
elseif mout == false && check == true
    [U,S,V] = deal([]);
elseif mout == true && check == false
    Strue = [];
elseif mout == true && check == true
end

