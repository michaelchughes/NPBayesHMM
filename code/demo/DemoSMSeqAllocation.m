% This demo shows a sequence of split or merge moves
%  with the full proposal + Met-Hastings accept ratio machinery in action
% DATA:  50 sequences of Gaussian observations with 2 dimensions
%           each 100 timesteps long
% PROCEDURE:
%  We start with config set to the ground truth via @initBPHMMCheat
%  User selects whether to perform SPLITS or MERGES via flag "doTrySplit"
%  For several trials, we repeatedly:
%     construct a "broken" input config at random (like in other demos)
%     apply the chosen move to this "broken" configuration
%     examine the accept-reject ratio that comes out
% OUTPUT
%   Figure that shows "BEFORE" and "AFTER" emission param configurations
%    for several trials of the split-merge MH proposal
%   Both accepts and rejects are shown.
% Note: we purposefully force the move to act on the "broken" behavior IDs,
%   otherwise this would be a less interesting demo

clear all;
doTrySplit = 1;
doPlot = 1;

Ktrue = 8;
if ~exist('data','var')
    data = genToySeqData_Gaussian( Ktrue, 2, 50, 100 );
end

algP.SM.featSelectDistr = 'splitBias+margLik';
algP.SM.doSeqUpdateThetaHatOnMerge = 0;

% Initialize to cheating configuration
Psi = initBPHMMCheat( data, defaultModelParams_BPHMM(data) );

tic;
for trial = 1:5
    brokePsi = Psi;
    if doTrySplit
        % ------------------------------- Construct broken "merged" config
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
        % ------------------------------- Construct broken "split" config
        origFeatIDs = randsample( 1:Ktrue, 1);
        kk = origFeatIDs(1);
        anchorIDs=[];
        for ii = 1:size(Psi.F,1)
            ttINDS = Psi.stateSeq(ii).z == origFeatIDs(1);
            if sum(ttINDS)>0
                Xseq = data.seq(ii);
                %Always assign to both
                brokePsi.F( ii, [kk Ktrue+1] ) = [1 1];
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
    % Sample the HMM trans + theta params for the broken config
    brokePsi.ThetaM = brokePsi.ThetaM.sampleAllTheta( data, brokePsi.stateSeq );
    brokePsi.TransM = brokePsi.TransM.sampleAllEta( brokePsi.F, brokePsi.stateSeq );
    [brokePsi, keepFeatIDs] = reallocateFeatIDs( brokePsi );
    
    % ------------------------------- Perform SM Met-Hastings move!
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


% NOTE: if instead we have more carefully crafted
%  split, for the brokePsi, our accept rates are miserable
% Because this kind of split is SUPER unlikely,
%  and the Hastings factor in accept ratio drowns out the gain in
%  likelihood
%                 aIDs = Xseq(1,ttINDS) < Psi.ThetaM.theta(kk).mu(1);
%                 bIDs = Xseq(1,ttINDS) >= Psi.ThetaM.theta(kk).mu(1);
%                 brokePsi.stateSeq(ii).z( ttINDS  ) = kk*aIDs + (Ktrue+1)*bIDs;
%                 if ~isempty(aIDs) && ii < 10
%                     anchorIDs(1) = ii;
%                 elseif ~isempty(bIDs)
%                     anchorIDs(2) = ii;
%                 end