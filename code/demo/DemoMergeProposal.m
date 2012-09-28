% This demo shows constructing a MERGE proposal
% DATA:  50 sequences of Gaussian observations with 2 dimensions
%           each 100 timesteps long
% PROCEDURE:
%  We start out with a configuration set to the ground truth
%     via @initBPHMMCheat
%  Next, we construct a "broken" input configuration by purposefully
%     splitting a single "true" behavior into two
%  Next, we apply the MERGE move to this "broken" configuration
%     to see if the proposal can come up with a better configuration
%     Note: we purposefully force merge at the particular "broken" behavior
%  Finally, we show that we can "reverse" this move
%     but constructing the exact broken input config via a Split
clear all;
close all;

% Initialize to cheating configuration
Ktrue = 4;
if ~exist('data','var')
    data = genToySeqData_Gaussian( Ktrue, 2, 50, 100 );
end
Psi = initBPHMMCheat( data, defaultModelParams_BPHMM(data) );

algP.SM.doSeqUpdateThetaHatOnMerge = 0;

% ======================================= Construct "broken" input
% Purposefully split feature 1,
%   and make sure the merge can recover these (more or less)
brokePsi = Psi;
for ii = 1:size(Psi.F,1)
    Xseq = data.seq(ii);
    ttINDS = Psi.stateSeq(ii).z==1;
    if ~isempty(ttINDS)
        aIDs = Xseq(1,ttINDS) < Psi.ThetaM.theta(1).mu(1);
        bIDs = Xseq(1,ttINDS) >= Psi.ThetaM.theta(1).mu(1);
        brokePsi.stateSeq(ii).z( ttINDS  ) = 1*aIDs + (Ktrue+1)*bIDs;
        brokePsi.F( ii, [1 Ktrue+1] ) = [~isempty(aIDs) ~isempty(bIDs)];
    end
end
brokePsi.ThetaM.K = Ktrue+1;
brokePsi.ThetaM = brokePsi.ThetaM.sampleAllTheta( data, brokePsi.stateSeq );
brokePsi.TransM = brokePsi.TransM.sampleAllEta( brokePsi.F, brokePsi.stateSeq );

% ======================================= Construct proposal!
% First, identify which feature to split (corresp. to feat 1 in orig Psi)
mergeFeats(1) = 1;
mergeFeats(2) = Ktrue+1;
candidates = find( Psi.F(:,1) >0 );

brokePsi.activeFeatIDs = mergeFeats;
anchorObjIDs = candidates( randperm( length(candidates) ) );
anchorObjIDs = anchorObjIDs(1:2);

% Actually construct proposal
starPsi = sampleMergeConfig( brokePsi, data, anchorObjIDs, brokePsi.activeFeatIDs, algP );

% Consider reverse move back to original (broken) config
[origPsi, logQ] = sampleSplitConfig( starPsi, data, anchorObjIDs, starPsi.activeFeatIDs, algP, brokePsi );

% ======================================= Plot input, proposal, & rev. move
figure(101);

subplot(1,3, 1);
plotEmissionParams( brokePsi, data );
title( 'BADLY SPLIT INPUT' );

starPsi = reallocateFeatIDs(starPsi);
subplot(1,3, 2);
plotEmissionParams( starPsi, data );
title( 'PROPOSED MERGE CONFIG' );

origPsi = reallocateFeatIDs(origPsi);
subplot(1,3, 3);
plotEmissionParams( origPsi, data );
title( 'REVERSE SPLIT BACK TO INPUT' );


% MAKE SURE WE ACTUALLY GOT BACK THE ORIGINAL CONFIG VIA REVERSE MOVE
Zall1 = horzcat(origPsi.stateSeq(:).z );
Zall2 = horzcat(brokePsi.stateSeq(:).z );
assert( all(Zall1==Zall2), 'ERROR!' );