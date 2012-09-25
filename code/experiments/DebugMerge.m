%clear all;

% Initialize to cheating configuration
Ktrue = 4;
if ~exist('data','var')
    data = genToySeqData_Gaussian( Ktrue, 2, 50, 100 );
end
Psi = initBPHMMCheat( data, defaultModelParams_BPHMM(data) );

%ChainHist = runBPHMM( {'SynthGaussian', 'nObj', 100}, {}, {1,1}, {'Niter', 100}, {'F.nTotal', 8} );
algP.SM.doSeqUpdateThetaHatOnMerge = 0;

% Purposefully split feature 1,
%   and make sure the merge can recover these (more or less)
brokePsi = Psi;
for ii = 1:size(Psi.F,1)
    if rand < 0.5   
        ttINDS = Psi.stateSeq(ii).z==1;
        brokePsi.stateSeq(ii).z( ttINDS  ) = Ktrue+1;
        brokePsi.F( ii, [1 Ktrue+1] ) = [0 1];
    end
end
brokePsi.ThetaM.K = Ktrue+1;
brokePsi.ThetaM = brokePsi.ThetaM.sampleAllTheta( data, brokePsi.stateSeq );
brokePsi.TransM = brokePsi.TransM.sampleAllEta( brokePsi.F, brokePsi.stateSeq );


% Perform a split move!
% First, identify which feature to split (corresp. to feat 1 in orig Psi)
mergeFeats(1) = 1;
mergeFeats(2) = Ktrue+1;
candidates = find( Psi.F(:,1) >0 );

brokePsi.activeFeatIDs = mergeFeats;
%[brokePsi, keepFeatIDs] = reallocateFeatIDs( brokePsi );

brokePsiIN = brokePsi;
tic;
for trial = 1:50
    anchorObjIDs = candidates( randperm( length(candidates) ) );
    anchorObjIDs = anchorObjIDs(1:2);
    
    %brokePsi = relabelFeatIDs( brokePsiIN, randperm(Ktrue+1) );
    
    starPsi = sampleMergeConfig( brokePsi, data, anchorObjIDs, brokePsi.activeFeatIDs, algP );
   
    % Consider reverse probability!
    [origPsi, logQ] = sampleSplitConfig( starPsi, data, anchorObjIDs, starPsi.activeFeatIDs, algP, brokePsi );
    % MAKE SURE WE ACTUALLY GET THE ORIGINAL CONFIG
    Zall1 = horzcat(origPsi.stateSeq(:).z );
    Zall2 = horzcat(brokePsi.stateSeq(:).z );
    assert( all(Zall1==Zall2), 'ERROR!' );
    
end
toc
return;
 figure(101);
 
    subplot(1,3, 1); 
    plotEmissionParams( brokePsi );
    
    starPsi = reallocateFeatIDs(starPsi);
    subplot(1,3, 2); 
    plotEmissionParams( starPsi );
    
    origPsi = reallocateFeatIDs(origPsi);
    subplot(1,3, 3);
    plotEmissionParams( origPsi );
    