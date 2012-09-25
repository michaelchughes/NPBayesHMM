%function [mu,sigma,a] = matrixNormalToNormal(M,V,K)
%
% Converts the parameters for a matrix normal A ~ MN(M,V,K) 
% into a  multivariate normal  A(:) ~ N(mu,sigma)
%
function [mu,sqrtsigma] = matrixNormalToNormal(M,sqrtV,sqrtinvK)

 mu = M(:);
 sqrtsigma = kron(sqrtinvK,sqrtV);
 
 % sigma = sqrtsigma'*sqrtsigma;
 
 % From Minka's paper, but order is wrong:
 %sigma = kron(V,inv(K));
 


