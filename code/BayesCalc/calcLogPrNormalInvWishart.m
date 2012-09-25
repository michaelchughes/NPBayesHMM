function logPr = calcLogPrNormalInvWishart( Mu, invSigma, PP )
% Mu, Sigma ~ NormInvWish( PP params )
% PP has fields
%  .precMu : precision of mean  (k in Murphy)
%  .Mu     :  mean of mean      (Mu in Murphy)
%  .ScaleMatrix :  DxD matrix   ( S  in Murphy)
%  .degFree :  degrees of freedom (v in Murphy)

logPrSigma = calcLogPrInvWishart( invSigma, PP );

logPrMu = calcLogPrGaussian( Mu, PP.mu, PP.precMu*invSigma );

logPr = logPrSigma + logPrMu;
end
