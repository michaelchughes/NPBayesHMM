function [splitPsi, logQ] = sampleSplitConfig( Psi, data, anchorObjIDs, featIDs, algParams, TargetPsi )
% Sample a proposal config of (F,z) that has one extra state than input
%INPUT:
%  Psi : input model config (contains F,z)
%  data : data object
%  anchorObjIDs : length 2 vector giving ids of anchor objects for proposal
%  featIDs : integer id of feature in Psi to split into two
%  TargetPsi : target model config. 
%    when provided, we don't actually sample, but instead calculate
%    probability "logQ" of sampling the target given the input Psi config
%OUTPUT:
%  splitPsi : output configuration (proposal)
%    has all its observed data sufficent statistics pre-calculated
%    for easy use in later calculations (accept ratios, etc.
%  logQ     : struct of transition probabilities (in log space)

% ------------------------------------------------------- INTERP INPUT
anchorObjIDs = anchorObjIDs(:)';
ii = anchorObjIDs(1);
jj = anchorObjIDs(2);

kold = featIDs(1);
otherObjIDs = setdiff( find( Psi.F( :, kold ) > 0   )', anchorObjIDs );
activeObjIDs = [anchorObjIDs otherObjIDs];

if exist( 'TargetPsi', 'var' ) && ~isempty( 'TargetPsi' )
   kA = TargetPsi.activeFeatIDs(1);
   kB = TargetPsi.activeFeatIDs(2);
else
    TargetPsi = [];
    kA = kold;
    kB = size(Psi.F,2) + 1;
end

% ------------------------------------------------------- INIT PROPOSAL
propF = Psi.F == 1;
propF(:, kold) = 0;
propF( ii, kA) = 1;
propF( jj, kB) = 1;
propFeatIDs = [kA kB];

propStateSeq = Psi.stateSeq;
for aa = 1:data.N
    propStateSeq( aa ).z( Psi.stateSeq(aa).z == kold ) = 0;
end
propStateSeq( ii ).z(  Psi.stateSeq( ii ).z == kold )  = kA;
propStateSeq( jj ).z(  Psi.stateSeq( jj ).z == kold )  = kB;

% Create deterministic trans params Eta
propTransM = Psi.TransM.getAllEta_PriorMean( activeObjIDs, propF, propFeatIDs );

% Create deterministic Theta for ALL features: others + [kA kB]
propThetaM = Psi.ThetaM.getAllTheta_PosteriorMean( data, propStateSeq, propFeatIDs, featIDs);
assert( propThetaM.K == length(propThetaM.Xstats), 'Badness');

% Init suff stats for F
featCounts = sum( propF(:,propFeatIDs) , 1 );
nObj = size(propF, 1) - length(otherObjIDs);
logQ.F = 0;
logQ.z = 0;

proposalsON = false(1, size(propF,2) );
proposalsON(propFeatIDs) = 1;
for aa = [ otherObjIDs( randperm(length(otherObjIDs)) ) ii jj]
    
    if aa == ii || aa == jj
        % Make sure to update F suff. stats for anchors
        %  since we've already "seen" them before, unlike other sequences
        featCounts = featCounts - propF(aa, propFeatIDs);
        nObj = nObj - 1;
        propThetaM = propThetaM.decXStats( aa, data, propStateSeq, propFeatIDs );
    end
    
    % ====================================================  Sample F(aa,:)
    % Calc Soft Evidence for the current sequence, including both new feats
    allActiveFeatIDs = find( propF(aa,:) | proposalsON );
    seqSoftEv = propThetaM.calcLogSoftEv( aa, data, allActiveFeatIDs );
    
    % Propose new feature assignments for the cur. seq. 
    priorPrFeats = featCounts ./ ( nObj + Psi.bpM.c );
    [newFaa, logQ_Faa] = sampleTwoSharedFeatsForSeq_Gibbs( aa, propFeatIDs, seqSoftEv, propTransM.seq(aa), priorPrFeats, anchorObjIDs, TargetPsi);
    propF(aa, propFeatIDs) = newFaa==1;
    logQ.F = logQ.F + logQ_Faa;    

    % ===============================================  Sample stateSeq(aa)
    seqTransM = propTransM.getSeqEtaWithFeats( aa, propF(aa,:) );
    propTransM = propTransM.setEta( aa, propF(aa,:), seqTransM.eta );
    [propStateSeq(aa).z, logQ_zaa] = sampleSingleStateSeq_WithSoftEv( aa, seqTransM, seqSoftEv(propF(aa,:), :), TargetPsi);
    logQ.z = logQ.z + logQ_zaa;
    
    % Update suff. stats for F
    featCounts = featCounts + propF(aa, propFeatIDs);
    nObj = nObj + 1;
    
    % Update Emission Params for new proposed features [kA, kB]
    %   using the newly allocated data from the current sequence
    propThetaM = propThetaM.incXStats( aa, data, propStateSeq, propFeatIDs );
    for kkID = 1:length(propFeatIDs)
        kk = propFeatIDs(kkID);
        PN = propThetaM.getPosteriorParams( propThetaM.Xstats(kk) );
        propThetaM.theta(kk) = propThetaM.getTheta_Mean( PN );
    end
end

if isempty( TargetPsi )
    propThetaM = propThetaM.updateAllXSuffStats( horzcat(propStateSeq(:).z), data );
end

splitPsi.F = propF;
splitPsi.stateSeq = propStateSeq;
splitPsi.ThetaM = propThetaM;
splitPsi.TransM = propTransM;
splitPsi.bpM = Psi.bpM;
splitPsi.activeFeatIDs = propFeatIDs;
splitPsi.anchorIDs  = anchorObjIDs;
