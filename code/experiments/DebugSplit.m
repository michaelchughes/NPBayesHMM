clear all;

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

brokePsiIN = brokePsi;
for trial = 1:50
    anchorObjIDs = candidates( randperm( length(candidates) ) );
    anchorObjIDs = anchorObjIDs(1:2);
    
    brokePsi = relabelFeatIDs( brokePsiIN, randperm(Ktrue-1) );
    
    [starPsi, logQ] = sampleSplitConfig( brokePsi, data, anchorObjIDs, brokePsi.activeFeatIDs, algP );
   
    % Consider reverse probability!
    [origPsi, logQREV] = sampleMergeConfig( starPsi, data, anchorObjIDs, starPsi.activeFeatIDs, algP, brokePsi );
    % MAKE SURE WE ACTUALLY GET THE ORIGINAL CONFIG
    Zall1 = horzcat(origPsi.stateSeq(:).z );
    Zall2 = horzcat(brokePsi.stateSeq(:).z );
    assert( all(Zall1==Zall2), 'ERROR!' );
end



figure(101);
 
subplot(1,3, 1); 
plotEmissionParams( brokePsi );
    
subplot(1,3, 2); 
plotEmissionParams( starPsi );
    
subplot(1, 3, 3);
plotEmissionParams( origPsi );