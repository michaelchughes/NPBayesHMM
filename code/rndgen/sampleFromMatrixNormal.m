%function S = sampleFromMatrixNormal(M,V,K,nSamples=1)
function S = sampleFromMatrixNormal(M,sqrtV,sqrtinvK,nSamples)

if ~exist('nSamples','var'), nSamples = 1; end

[mu,sqrtsigma] = matrixNormalToNormal(M,sqrtV,sqrtinvK);

S = mu + sqrtsigma'*randn(length(mu),1);
S = reshape(S,size(M));
