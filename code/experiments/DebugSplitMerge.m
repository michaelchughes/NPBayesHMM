clear all;
doTrySplit = 1;
doPlot = 0;

Ktrue = 8;
if ~exist('data','var')
    data = genToySeqData_Gaussian( Ktrue, 2, 50, 100 );
end

algP.SM.featSelectDistr = 'splitBias+margLik';
algP.SM.doSeqUpdateThetaHatOnMerge = 0;

% Initialize to cheating configuration
Psi = initBPHMMCheat( data, defaultModelParams_BPHMM(data) );

% Construct a configuration where two "true" features are merged into one
tic;
for trial = 1:10
    brokePsi = Psi;
    if doTrySplit
        origFeatIDs = randsample( 1:Ktrue, 2);
        
        for ii = 1:size(Psi.F,1)
            ttINDS = Psi.stateSeq(ii).z == origFeatIDs(1) | Psi.stateSeq(ii).z == origFeatIDs(2);
            if sum(ttINDS)>0
                brokePsi.F( ii, origFeatIDs ) = [1 0];
                brokePsi.stateSeq(ii).z( ttINDS ) = origFeatIDs(1);
            end
        end
        brokePsi.ThetaM.K = Ktrue-1;
        
        brokePsi.activeFeatIDs = origFeatIDs(1);
        brokePsi.anchorIDs = find( brokePsi.F(:,origFeatIDs(1) )>0, 2);
    else        
        origFeatIDs = randsample( 1:Ktrue, 1);
        
        anchorIDs=[];
        for ii = 1:size(Psi.F,1)
            ttINDS = Psi.stateSeq(ii).z == origFeatIDs(1);
            if sum(ttINDS)>0
                brokePsi.F( ii, [origFeatIDs(1) Ktrue+1] ) = [1 1];
                if rand < 0.5
                    anchorIDs(1)=ii;
                    brokePsi.stateSeq(ii).z( ttINDS ) = Ktrue+1; 
                else
                    anchorIDs(2)=ii;
                end
            end
        end
        assert( all(anchorIDs~=0) && length(anchorIDs)==2,'BADNESS');
        brokePsi.ThetaM.K = Ktrue+1;
        
        brokePsi.activeFeatIDs = [origFeatIDs(1) Ktrue+1];
        brokePsi.anchorIDs = anchorIDs;
    end
    
    brokePsi.ThetaM = brokePsi.ThetaM.sampleAllTheta( data, brokePsi.stateSeq );
    brokePsi.TransM = brokePsi.TransM.sampleAllEta( brokePsi.F, brokePsi.stateSeq );
    
    [brokePsi, keepFeatIDs] = reallocateFeatIDs( brokePsi );

    
[starPsi, Stats] = sampleSplitMerge_SeqAlloc( brokePsi, data, algP );
if Stats.nAccept==1
    outcomeStr = 'accept';
else
    outcomeStr = 'reject';
end

fprintf( ' %s : rho=%.2f\n', outcomeStr, Stats.rho );

if doPlot
figure;
subplot(1,2,1);
plotEmissionParams( brokePsi, data);
title('BEFORE');

subplot(1,2,2);
plotEmissionParams( starPsi, data);
title(['AFTER: ' outcomeStr]);
end
end
fprintf( '............. elapsed time %.1f\n', toc );