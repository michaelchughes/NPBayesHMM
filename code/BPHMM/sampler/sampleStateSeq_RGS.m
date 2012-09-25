function [z, logQ] = sampleStateSeq_RGS( logSoftEv, eta_ii, F_ii, ii, TargetPsi )
% Draw new assignments for the discrete hidden state sequence
%  for each time series object, under *restricted settings*.
% Uses a fast message backward, sample forwards algorithm.

ks = find( F_ii > 0 );
T  = size(logSoftEv,2);
K  = length(ks);

if K == 1
    z = ks(1)*ones( 1, Tii );
    continue;
end

pi_init = 1/K*ones(1,K);
pi      = eta_ii;

normC = max( logSoftEv, [], 1);
logSoftEv = bsxfun( @minus, logSoftEv, normC );
Lik = exp( logSoftEv );
     
SEED = randomseed();
randomseed( SEED+1 );
[z, logQ] = SampleHMMStateSeqWithQsC( pi, Lik, pi_init, SEED(1) );
z = ks(z);
   
end % main function


