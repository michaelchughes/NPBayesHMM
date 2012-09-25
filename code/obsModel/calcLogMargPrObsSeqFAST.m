function margLogPr = calcLogMargPrObsSeqFAST( LL, Eta )
% calcLogMargPrObsSeqFAST
% Provides fast calculation of marginal likelihood
%    for a particular sequence of data, given HMM parameters in the form
%      soft evidence LL, where LL(kk,tt) = p( X{ii}(tt) | theta(kk) )
%      transition weights Pz, which may be non-normalized
% Massages input data and then calls a super-efficient MEX function
%   "FilterFwdC" to perform dynamic programming.
%INPUT
% LL  := KiixTii matrix of log likelihoods soft evidence
% Pz  := KiixKii matrix of transition weights (non-normalized)
%OUTPUT
%  margLogPr  := scalar value of log( X_ii | F_ii, theta, Eta_ii ) 

K = size(Eta,2);
if K == 0
    margLogPr = -Inf;
    return;
end

Pi = bsxfun( @rdivide, Eta, sum(Eta,2) );

% Need to turn log_lik into lik.  We know lik = exp( log_lik )
%   For numerical stability, we find M = max( log_lik )
%       we compute L = exp( log_lik - M  )
%       and we thus have lik up to multiplicative constant
%              since lik \propto L = exp( log_lik ) / exp( M )
%   We find a unique "normalizer" M_t  for each time step
   
% normC = 1 x T
normC = max( LL,[],1);
    
if K == 1
    margLogPr = sum(normC );
    return;
end
Lik = exp( bsxfun( @minus, LL, normC ) );
    
[~,margLogPr] = FilterFwdC( Pi, Lik, 1/K*ones(1,K) );
margLogPr = margLogPr + sum( normC );

end % main function

