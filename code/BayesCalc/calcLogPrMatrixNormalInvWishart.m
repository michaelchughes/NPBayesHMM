function [logPr] = calcLogPrMatrixNormalInvWishart( A, invSigma, PP )

logPrSigma = calcLogPrInvWishart( invSigma, PP );

cholInvSigma = chol( invSigma );
logDetInvSigma = 2*sum( log( diag( cholInvSigma ) ) );
logDetSigma    = -logDetInvSigma; % det(S) = 1/det(inv(S))

v = PP.degFree;
SM  = PP.ScaleMat;
CC = PP.invAScaleMat;
[D DR] = size( CC );

if isfield( PP, 'MeanMat') && ~isempty( PP.MeanMat )
    UU = A - PP.MeanMat;
else
    UU = A; % assume zero mean
end

logPrA = 0.5*D*log( det(CC) ) - 0.5*DR*D*log( 2*pi ) ...
          - 0.5*DR*logDetSigma ...
          - 0.5*trace( UU'*invSigma*UU*CC );
      
      
logPr = logPrA + logPrSigma;
