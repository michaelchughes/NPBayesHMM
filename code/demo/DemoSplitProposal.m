% This demo shows constructing a split proposal
% DATA:  50 sequences of Gaussian observations with 2 dimensions
%           each 100 timesteps long
% PROCEDURE:
%  We start out with a configuration set to the ground truth
%     via @initBPHMMCheat
%  Next, we construct a "broken" input configuration by purposefully
%     merging two ground truth sequences into one
%  Next, we apply the split move to this "broken" configuration
%     to see if the proposal makes sense
%     Note: we purposefully force split at the particular "broken" behavior
%  Finally, we show that we can "reverse" this move
%     but constructing the exact broken input config via a Merge

clear all;
close all;

Ktrue = 4;
if ~exist('data','var')
    data = genToySeqData_Gaussian( Ktrue, 2, 50, 100 );
end

algP.SM.doSeqUpdateThetaHatOnMerge = 0;

% Initialize to cheating configuration
Psi = initBPHMMCheat( data, defaultModelParams_BPHMM(data) );

% Purposefully merge features 1 and 2,
%   and make sure the split move can recover these
brokePsi = Psi;
for ii = 1:size(Psi.F,1)
    ttINDS = Psi.stateSeq(ii).z == 1 | Psi.stateSeq(ii).z == 2;
    if sum(ttINDS)>0
        brokePsi.F( ii, [1 2] ) = [1 0];
        brokePsi.stateSeq(ii).z( ttINDS ) = 1;
    end
end
brokePsi.ThetaM.K = Ktrue+1;
brokePsi.ThetaM = brokePsi.ThetaM.sampleAllTheta( data, brokePsi.stateSeq );
brokePsi.TransM = brokePsi.TransM.sampleAllEta( brokePsi.F, brokePsi.stateSeq );

% Perform a split move!
% First, identify which feature to split (corresp. to feat 1 in orig Psi)
k2split = 1;

candidates = find( brokePsi.F(:,1) >0 );
brokePsi.activeFeatIDs = k2split;
[brokePsi, keepFeatIDs] = reallocateFeatIDs( brokePsi );

anchorObjIDs = candidates( randperm( length(candidates) ) );
anchorObjIDs = anchorObjIDs(1:2);


[propPsi, logQ] = sampleSplitConfig( brokePsi, data, anchorObjIDs, brokePsi.activeFeatIDs, algP );

% Consider reverse probability!
[origPsi, logQREV] = sampleMergeConfig( propPsi, data, anchorObjIDs, propPsi.activeFeatIDs, algP, brokePsi );
% MAKE SURE WE ACTUALLY GET THE ORIGINAL CONFIG
Zall1 = horzcat(origPsi.stateSeq(:).z );
Zall2 = horzcat(brokePsi.stateSeq(:).z );
assert( all(Zall1==Zall2), 'ERROR!' );

figure(101);

subplot(1,3, 1);
plotEmissionParams( brokePsi, data );
title( 'BADLY MERGED INPUT' );

subplot(1,3, 2);
plotEmissionParams( propPsi, data );
title( 'PROPOSED SPLIT CONFIG' );

subplot(1, 3, 3);
% Reallocate to remove the superfluous feature left over after the merge
origPsi = reallocateFeatIDs( origPsi );
plotEmissionParams( origPsi, data );
title( 'REVERSE MERGE BACK TO INPUT' );