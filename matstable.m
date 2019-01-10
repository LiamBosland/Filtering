function [stable] = matstable(A,ploteig)
% This functions checks if all eigenvalues of the provided matrix A are
% within the unit circle (thus matrix A is stable) and, if provided by
% LOGICAL 'ploteig', plots the eigenvalues versus the unit circle.

if nargin == 1 % Only matrix A is provided as input, eigenvalues are not plotted
    ploteig = false;
end

E = eig(A); r = real(E); im = imag(E);
s = 0;
for k = 1 : length(E)
    if abs(r(k)) <= 1 && abs(im(k)) <= 1
        s = s + 1;
    end
end
if s == length(E)
    stable = true;
elseif s ~= length(E)
    stable = false;
end
if ploteig == true
    [xc,yc] = deal(linspace(0,2*pi,length(E)));
    figure;
    plot(sin(xc),cos(yc),'k-',im,r,'ro');
    xlabel('Imaginary axis'); ylabel('Real axis'); title('Eigenvalues of matrix A*');
    legend('unit circle','eigenvalues A*'); grid on;
end
end

