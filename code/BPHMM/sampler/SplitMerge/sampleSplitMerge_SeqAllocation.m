function [F,stateSeq, Stats] = sampleSplitMerge_SeqAllocation( F, stateSeq, data_struct, hyperparams, model, settings, Stats, Debug )
% Split Merge with Sequential Allocation updates

% Create structure for capturing descriptive stats about this move
if ~isfield( Stats, 'nAccept' )
    M = settings.F.nSMTrials;
    Stats.nAccept = 0;
    Stats.nTrial  = 0;
    Stats.avgRatio = 0;
    
    Stats.rhos = zeros( 1, M );
    Stats.moveTypes = zeros( 1, M );
    Stats.moveSizes = zeros( 1, M );
    Stats.anchorObjIDs = zeros( M, 2);
    Stats.featIDs = zeros( M, 2);
    
    Stats.SPLIT.nAccept = 0;
    Stats.SPLIT.nTrial = 0;
    Stats.MERGE.nAccept = 0;
    Stats.MERGE.nTrial  = 0;
end


Xstats_CUR = getXSuffStats( F, stateSeq, data_struct, model );
if exist( 'Debug', 'var' ) && isfield( Debug, 'anchorObjIDs' )
    anchorObjIDs = Debug.anchorObjIDs;
    ii = Debug.anchorObjIDs(1);
    jj = Debug.anchorObjIDs(2);
    
    featIDs = Debug.featIDs;
    ki = Debug.featIDs(1);
    kj = Debug.featIDs(2);
    assert( F(ii,ki) == 1, 'BAD CHOICE');
    assert( F(jj,kj) == 1, 'BAD CHOICE');
    ps = getSplitMergeFeatChoiceProbsBetter( F, stateSeq, data_struct, model, jj, ki, Xstats_CUR );
else
    % Randomly sample anchor objects
    anchorObjIDs = randsample( 1:size(F,1), 2 );
    ii = anchorObjIDs(1);
    jj = anchorObjIDs(2);

    ki = multinomial_single_draw( F(ii,:) );
    %ps = getSplitMergeFeatChoiceProbs( F, jj, ki );
    %ps = getSplitMergeFeatChoiceProbs( F, stateSeq, data_struct, model, jj, ki );
    ps = getSplitMergeFeatChoiceProbsBetter( F, stateSeq, data_struct, model, jj, ki, Xstats_CUR );
    kj = multinomial_single_draw( ps );
    
    assert( ~isempty( kj ), 'Crap: Empty kj' );
    featIDs = [ki kj];
end

Psi.F = F;
Psi.stateSeq = stateSeq;

if ki == kj
    % =========================================== SPLIT
    moveDescrStr = 'SPLIT';

    [propPsi, logQ] = getSplitConfig_SeqAllocation( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs );

    origPsi = Psi;
    origPsi.activeFeatIDs = ki;
    
    [origPsi, logQ_Rev] = getMergeConfig_SeqAllocation( propPsi, data_struct, hyperparams, model, anchorObjIDs, propPsi.activeFeatIDs, origPsi );
else
    % =========================================== MERGE
    moveDescrStr = 'MERGE';
    origPsi = Psi;
    origPsi.activeFeatIDs = [ki kj];

    [propPsi, logQ] = getMergeConfig_SeqAllocation( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs );
    
    [origPsi, logQ_Rev] = getSplitConfig_SeqAllocation( propPsi, data_struct, hyperparams, model, anchorObjIDs, propPsi.activeFeatIDs, origPsi );    
    
end

% Compute probability of choosing ki, kj for FORWARD MOVE
qFeats = ps; %getSplitMergeFeatChoiceProbs( Psi.F, Psi.stateSeq, data_struct, model, jj, ki );
logQ.featChoice = log( 1/sum( Psi.F(anchorObjIDs(1),:) ) ) ...
    + log( qFeats(kj)/sum(qFeats)  );

% Compute probability of choosing corresponding ki, kj for REVERSE MOVE
Xstats_PROP = getXSuffStats( propPsi.F, propPsi.stateSeq, data_struct, model );

if length(  propPsi.activeFeatIDs ) == 2
    propki = propPsi.activeFeatIDs(1);
    propkj = propPsi.activeFeatIDs(2);
else
    propki = propPsi.activeFeatIDs(1);
    propkj = propPsi.activeFeatIDs(1);
end
%qRevFeats = getSplitMergeFeatChoiceProbs(  propPsi.F, propPsi.stateSeq, data_struct, model, jj, propki );
qRevFeats = getSplitMergeFeatChoiceProbsBetter(  propPsi.F, propPsi.stateSeq, data_struct, model, jj, propki, Xstats_PROP);
logQ_Rev.featChoice = log( 1/sum( propPsi.F(anchorObjIDs(1),:) ) ) ...
    + log( qRevFeats(propkj)/sum(qRevFeats)  );

% Total up probabilities of FORWARD (Q) and REVERSE (Q_Rev) moves
logQ.all = logQ.featChoice + logQ.F + logQ.z;
logQ_Rev.all = logQ_Rev.featChoice + logQ_Rev.F + logQ_Rev.z;

objIDs = find( sum( propPsi.F(:,propPsi.activeFeatIDs), 2 ) > 0 )';
logPr_Cur = calcLogJointPr_RestrictedSet_BPHMM( Psi, data_struct, hyperparams, model, objIDs, Xstats_CUR );
logPr_Prop = calcLogJointPr_RestrictedSet_BPHMM( propPsi, data_struct, hyperparams, model, objIDs, Xstats_PROP );


logPrAccept = logPr_Prop.all - logPr_Cur.all + logQ_Rev.all - logQ.all;
rho =  exp( logPrAccept );
rho = min(1, rho);
doAccept = rand < rho;

%ShowDebugInfoForSeqSplitMerge
if (  doAccept )
    newPsi = propPsi;
    Stats.nAccept = Stats.nAccept + 1;
else
    newPsi = Psi;
end

% Update descriptive stats with success/failure info from this move
rho = min(1, exp( logPrAccept ) );
Stats.avgRatio = ( Stats.nTrial*Stats.avgRatio + rho) / (Stats.nTrial+1);
Stats.nTrial = Stats.nTrial + 1;
Stats.rhos( Stats.nTrial ) = rho;
Stats.anchorObjIDs( Stats.nTrial, : ) = anchorObjIDs;
Stats.featIDs( Stats.nTrial, : ) = [ki kj];
Stats.moveSizes( Stats.nTrial )   = length( objIDs );
switch moveDescrStr
    case 'SPLIT'
        Stats.moveTypes( Stats.nTrial ) = 0;
    case 'MERGE'
        Stats.moveTypes( Stats.nTrial ) = 1;
end
Stats.( moveDescrStr ).nAccept = Stats.( moveDescrStr ).nAccept + doAccept;
Stats.( moveDescrStr ).nTrial = Stats.( moveDescrStr ).nTrial + 1;

% Remove empty columns of F, and rename state sequence appropriately
colSums = sum( newPsi.F, 1);
keepIDs = find( colSums > 0 );
stateSeq = newPsi.stateSeq;
for ii = 1:length( stateSeq )
    for kk = 1:length( keepIDs )
        stateSeq(ii).z( stateSeq(ii).z == keepIDs(kk)  ) = kk;
    end
end
F = newPsi.F( :, keepIDs );

end % main function
