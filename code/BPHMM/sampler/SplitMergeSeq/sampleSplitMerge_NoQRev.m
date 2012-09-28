function [newPsi, Stats] = sampleSplitMerge_NoQRev( Psi, data, algParams )
% Split Merge via Not-Quite-Valid Sequential Allocation
% Proposes new candidate configuration for feat matrix F and stateSeq z
%  and accepts or rejects via Metropolis-Hastings
% HOWEVER, this *does not* account for the probability of the reverse move
%  in the acceptance ratio... there is still an acceptance ratio, but
%  it is entirely based on the joint probability under the model,
%  not on the proposal probabilities.
% Proposal selects anchor objects and features to split/merge at random,
%   unless provided (see DEBUG NOTE below)
%INPUT:
%  Psi  : model config input
%  data : SeqData object
%  algParams : relevant params for this method
%    ** SM.featSelectDistr can be one of          
%       'random' : choose features to merge at random
%       'splitBias+margLik' (default) : use smart criteria to
%                    make successful proposals
%    ** SM.doSeqUpdateThetaHatOnMerge
%       if 1, refine Theta after visiting every seq. during merge proposal
%          0 (default), just estimate Theta initially and stick with it
%OUTPUT:
%  newPsi : model configuration as result of MH proposal
%    equals a new configuration when proposal accepted
%    or the input config Psi when rejected
%  Stats  : struct indicating type of move and whether accepted/rejected
%ASSUMES:
%  Psi comes with correct sufficient statistics for observed data (Xstats)
%DEBUG NOTE:
%   Can provide input "Psi" with fields "anchorIDs" and "activeFeatIDs"
%     to indicatea specific proposal to try. 
%   E.g. to merge features 1 and 5 using objects 101 and 202 as seeds, set
%    Psi.anchorIDs = [101 202];
%    Psi.activeFeatIDs = [1 5];  Note F(101,1) and F(202,5) must be ON.

if isfield( Psi, 'activeFeatIDs' )
    [anchorIDs, featIDs, qFWD] = sampleAnchorAndFeatIDsToSplitMerge( Psi, data, algParams );
else
    [anchorIDs, featIDs, qFWD] = sampleAnchorAndFeatIDsToSplitMerge( Psi, data, algParams );
    Psi.anchorIDs        = anchorIDs;
    Psi.activeFeatIDs = featIDs;
end

if featIDs(1) == featIDs(2)
    % =========================================== SPLIT
    moveDescrStr = 'ADD';
    [propPsi] = sampleSplitConfig( Psi, data, anchorIDs, featIDs(1), algParams );
    %[propPsi, logQ] = sampleSplitConfig( Psi, data, anchorIDs, featIDs(1), algParams );
    %[~, logQ_Rev]   = sampleMergeConfig( propPsi, data, anchorIDs, propPsi.activeFeatIDs, algParams, Psi );
else
    % =========================================== MERGE
    moveDescrStr = 'DEL';   
    [propPsi] = sampleMergeConfig( Psi, data, anchorIDs, featIDs, algParams );
    
    %[propPsi, logQ] = sampleMergeConfig( Psi, data, anchorIDs, featIDs, algParams );
    %[~, logQ_Rev]   = sampleSplitConfig( propPsi, data, anchorIDs, propPsi.activeFeatIDs, algParams, Psi );    
end

% Total up probabilities of FORWARD (Q) and REVERSE (Q_Rev) moves
[~, ~, qREV] = sampleAnchorAndFeatIDsToSplitMerge( propPsi, data, algParams );
logQ_Rev.all = log(qREV); % + logQ_Rev.F + logQ_Rev.z;
logQ.all     = log(qFWD); % + logQ.F + logQ.z;

% NB: not passing data as an arg here
%   means that we trust the stored X suff stats in Psi! Yay efficiency.
logPr_Cur  = calcJointLogPr_BPHMMState( Psi );
logPr_Prop = calcJointLogPr_BPHMMState( propPsi );

logPrAccept = logPr_Prop.all - logPr_Cur.all + logQ_Rev.all - logQ.all;
rho = exp( logPrAccept );
rho = min(1, rho);
doAccept = rand < rho;

if (  doAccept )
    newPsi = propPsi;
    % Remove empty columns of F, and rename state sequence appropriately
    newPsi = reallocateFeatIDs( newPsi );
else
    newPsi = Psi;
end

% Strip off info used internally to identify how to construct proposals
newPsi = rmfield( newPsi, 'anchorIDs');
newPsi = rmfield( newPsi, 'activeFeatIDs');

Stats.nAccept = doAccept;
Stats.rho = rho;
Stats.moveDescr = moveDescrStr;

end % main function