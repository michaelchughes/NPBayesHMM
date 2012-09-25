function y = logMvGamma( x, D)
% Compute logarithm multivariate Gamma function.
%   operates on each entry of input matrix/vector x
%INPUT
%   x : any scalar/vector/matrix
%   D : integer dimension of the mv Gamma function
%MATH DETAILS -------------------------------------------------------
% Gamma_D(x) = pi^(D(D-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
% log Gamma_D(x) = D(D-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
% Credit: Michael Chen (sth4nth@gmail.com).

s = size(x);
ds = (1-(1:D)')/2;

% Force input x to be a row vector
X = reshape(x,1,prod(s));

% X(dd, dim) = input(dim) - ds(dd)
X = bsxfun(@plus, X, ds   );

y = D*(D-1)/4*log(pi) + sum(gammaln(X),1);

% Force output back to size of original input
y = reshape(y,s);
