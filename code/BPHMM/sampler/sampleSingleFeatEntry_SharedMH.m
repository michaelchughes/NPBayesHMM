function [newFik, logMargPrObs_New] = ... 
    sampleSingleFeatEntry_SharedMH( curFik, Mk, nObj, c0, propEta, propLL, logMargPrObs_Cur)
% Sample a binary indicator F(ii,kk), which determines whether
%         sequence/object ii possesses feature/behavior kk
% Uses metropolis-hastings proposals, likelihoods calc'd by forward-bkward.
%INPUT
%  curFik : current binary value of entry being sampled  
%  Mk     : # of other sequences possessing feature kk
%  nObj   : total # of sequences in the dataset
%  c0     : BP conc parameter, influences prior prob. of sharing features
%  propPz : proposed HMM transition params eta for current sequence
%             Kii x Kii matrix 
%  propLL : proposed soft evidence matrix 
%             Kii x T matrix
%  logMargPrObs_Cur  : p( x(ii) | F(ii,:), theta, eta ) for current config
%           the proposal considered here is compared to this, and 
%           accepted/rejected accordingly
%OUTPUT
%  newFik : new binary value of F(ii,kk)
%             = ~curFik if accepted
%             = curFik if rejected
%  logMargPrObs_New  : new p( x(ii) | F(ii,:), theta, eta ) value

% -------------------------------- Handle edge case
% Sometimes, we cannot sample an entry in F with this method
% This occurs if *only* the current seq. possesses this feature,
%   in which case we just need to return gracefully
if Mk == 0 || isempty( propEta )
    newFik = curFik;
    logMargPrObs_New = logMargPrObs_Cur;
    return;
end


% -------------------------------- Prior terms
%  Pr( F(ii,kk) = 1 | F(~ii,kk) ) \propto  m^{not i}_k / (c0 + N -1 )
if curFik == 1
    PrF_diff  = (nObj + c0 - 1 - Mk)/Mk ;
elseif curFik == 0
    PrF_diff  = Mk/(nObj + c0 - 1 - Mk) ;
end

% -------------------------------- Likelihood terms
logMargPrObs_Prop = calcLogMargPrObsSeqFAST( propLL, propEta );

% -------------------------------- Accept or Reject!
PrAccept = PrF_diff * exp( logMargPrObs_Prop - logMargPrObs_Cur );
if rand < PrAccept
    newFik = ~curFik;
    logMargPrObs_New = logMargPrObs_Prop;
else
    newFik = curFik;
    logMargPrObs_New = logMargPrObs_Cur;
end
