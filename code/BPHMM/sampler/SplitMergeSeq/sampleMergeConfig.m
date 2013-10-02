function [mergePsi, logQ] = sampleMergeConfig( Psi, data, anchorObjIDs, featIDs, algParams, TargetPsi )
% Sample a proposed config of (F,z) that merges two states into one
%INPUT:
%  Psi : input model config (contains F,z)
%  data : data object
%  anchorObjIDs : length 2 vector giving ids of anchor objects for proposal
%  featIDs : integer ids of features in Psi to merge into one
%  TargetPsi : target model config. 
%    when provided, we don't actually sample, but instead calculate
%    probability "logQ" of sampling the target from the input Psi config
%OUTPUT:
%  splitPsi : output configuration (proposal)
%    has all its observed data sufficent statistics pre-calculated
%    for easy use in later calculations (accept ratios, etc.
%  logQ     : struct of transition probabilities (in log space)



logQ.F   = 0;
logQ.z   = 0;

anchorObjIDs = anchorObjIDs(:)';
ii = anchorObjIDs(1);
jj = anchorObjIDs(2);

allObjIDs = find( sum( Psi.F(:,featIDs), 2 ) > 0 )';
otherObjIDs = setdiff( allObjIDs, anchorObjIDs );


if exist( 'TargetPsi', 'var' ) && ~isempty( 'TargetPsi' )
   kmerge = TargetPsi.activeFeatIDs(1);
else
    TargetPsi = [];
    kmerge = featIDs(1);
end

propF = Psi.F==1;
propF(:, featIDs(2) ) = 0;

propFeatIDs = kmerge;
propF(allObjIDs, kmerge) = 1;

propStateSeq = Psi.stateSeq;
if algParams.SM.doSeqUpdateThetaHatOnMerge
    for aa = 1:data.N
        z = Psi.stateSeq(aa).z;
        if aa == ii || aa == jj
            propStateSeq(aa).z( z==featIDs(1)|z==featIDs(2) ) = kmerge;
        else
            propStateSeq(aa).z( z==featIDs(1)|z==featIDs(2) ) = 0;
        end
    end
else
    for aa = 1:data.N
        aIDs =  propStateSeq(aa).z == featIDs(1);
        bIDs =  propStateSeq(aa).z == featIDs(2);
        propStateSeq(aa).z( aIDs | bIDs ) = kmerge;
    end
end

% Create deterministic trans params Eta
propTransM = Psi.TransM.getAllEta_PriorMean( allObjIDs, propF, propFeatIDs );

% Create deterministic Theta for ALL features: others + [kA kB]
propThetaM = Psi.ThetaM.getAllTheta_PosteriorMean( data, propStateSeq, propFeatIDs, featIDs );
assert( propThetaM.K == length(propThetaM.Xstats), 'Badness');

if ~exist( 'TargetPsi', 'var' )
    TargetPsi = [];
elseif exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
    TargetPsi.externalFeatIDs = kmerge;
end

proposalsON = false(1, size(propF,2) );
proposalsON(propFeatIDs) = 1;
for aa = [ otherObjIDs( randperm( length(otherObjIDs) ) ) anchorObjIDs]

    % Calc Soft Evidence for the current sequence, including both new feats
    allActiveFeatIDs = find( propF(aa,:) | proposalsON );
    seqSoftEv = propThetaM.calcLogSoftEv( aa, data, allActiveFeatIDs );
    
    
    % ==============================================  Sample stateSeq(aa)
    seqTransM = propTransM.getSeqEtaWithFeats( aa, propF(aa,:) );
    [propStateSeq(aa).z, logQ_zaa] = sampleSingleStateSeq_WithSoftEv( aa, seqTransM, seqSoftEv(propF(aa,:), :), TargetPsi);
    logQ.z = logQ.z + logQ_zaa;
    
    if algParams.SM.doSeqUpdateThetaHatOnMerge
        if aa == ii || aa == jj
           propThetaM = propThetaM.decXStats( aa, data, propStateSeq, propFeatIDs );
        end
        % Update Emission Params for new proposed features
        %   using the newly allocated data from the current sequence
        propThetaM = propThetaM.incXStats( aa, data, propStateSeq, propFeatIDs );
        for jj = 1:length(propFeatIDs)
            kk = propFeatIDs(jj);
            PN = propThetaM.getPosteriorParams( propThetaM.Xstats(kk) );
            propThetaM.theta(kk) = propThetaM.getTheta_Mean( PN );
        end
    end
end

if isempty(TargetPsi)
    propThetaM = propThetaM.updateAllXSuffStats( horzcat(propStateSeq(:).z), data );
end

mergePsi.F = propF;
mergePsi.stateSeq = propStateSeq;
mergePsi.ThetaM   = propThetaM;
mergePsi.TransM   = propTransM;
mergePsi.bpM = Psi.bpM;
mergePsi.activeFeatIDs = propFeatIDs;
mergePsi.anchorIDs  = anchorObjIDs;