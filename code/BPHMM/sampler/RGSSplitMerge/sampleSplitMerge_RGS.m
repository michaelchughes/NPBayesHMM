function [F,stateSeq,TS,theta, Stats] = sampleSplitMerge_RGS( F, stateSeq, TS, theta, data_struct, hyperparams, model, settings, Stats, Debug )
% Split Merge with Restricted Gibbs Scan updates
% This proposal instantiates *all* variables
%   feature matrix, state sequence, transition params, emission params


% Create structure for capturing descriptive stats about this move
if ~isfield( Stats, 'nAccept' )
    M = settings.F.nSMTrials;
    Stats.nAccept = 0;
    Stats.nTrial  = 0;
    
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
Psi.theta = theta;
Psi.TS = TS;

if ki == kj
    % =========================================== SPLIT
    moveDescrStr = 'SPLIT';
    Psi.activeFeatIDs = [ki kj];
    
    [launchPsi] = launchSplitMove_RGS( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs );
    for rr = 1:settings.F.nSweepsRGS;
        launchPsi = RestrictedGibbsScan_Split( launchPsi, data_struct, hyperparams, model, anchorObjIDs );
    end
    [propPsi, logQ] = RestrictedGibbsScan_Split( launchPsi, data_struct, hyperparams, model, anchorObjIDs );
    
    % Comple probability of reverse move:  propPsi >>>merge>>> Psi 
    mergePsi = launchMergeMove_RGS( propPsi, data_struct, hyperparams, model, propPsi.activeFeatIDs );
    for rr = 1: settings.F.nSweepsRGS;
        mergePsi = RestrictedGibbsScan_Merge( mergePsi, data_struct, hyperparams, model );
    end       
    [~, logQ_Rev] = RestrictedGibbsScan_Merge( mergePsi, data_struct, hyperparams, model, Psi );
    
else
    % =========================================== MERGE
    moveDescrStr = 'MERGE';
    Psi.activeFeatIDs = [ki kj];

    
    launchPsi = launchMergeMove_RGS( Psi, data_struct, hyperparams, model, featIDs );
    for rr = 1: settings.F.nSweepsRGS;
        launchPsi = RestrictedGibbsScan_Merge( launchPsi, data_struct, hyperparams, model );
    end       
    [propPsi, logQ] = RestrictedGibbsScan_Merge( launchPsi, data_struct, hyperparams, model );
    
    % Comple probability of reverse move:  propPsi >>>merge>>> Psi 
    [splitPsi] = launchSplitMove_RGS( propPsi, data_struct, hyperparams, model, anchorObjIDs, propPsi.activeFeatIDs );
    for rr = 1: settings.F.nSweepsRGS;
        splitPsi = RestrictedGibbsScan_Split( splitPsi, data_struct, hyperparams, model, anchorObjIDs );
    end
    [~, logQ_Rev] = RestrictedGibbsScan_Split( splitPsi, data_struct, hyperparams, model, anchorObjIDs, Psi );
    
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
logQ.all = logQ.featChoice + logQ.all;
logQ_Rev.all = logQ_Rev.featChoice + logQ_Rev.all;

objIDs = find( sum( propPsi.F(:,propPsi.activeFeatIDs), 2 ) > 0 )';
logPr_Cur = calcLogJointPr_RGS_BPHMM( Psi, data_struct, hyperparams, model, objIDs, Xstats_CUR );
logPr_Prop = calcLogJointPr_RGS_BPHMM( propPsi, data_struct, hyperparams, model, objIDs, Xstats_PROP );


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

if doAccept 
    propLL=calcLogJointPr_BPHMM( propPsi.F, propPsi.stateSeq, hyperparams, data_struct, model );
    curLL=calcLogJointPr_BPHMM( Psi.F, Psi.stateSeq, hyperparams, data_struct, model );
    %assert( curLL.all - propLL.all < 500, 'Bad Merger!' );
    if ( curLL.all - propLL.all ) > 500
        fprintf( 'WARNING: %s move caused big drop in collapsed log probability\n', moveDescrStr );
    end

end

% Update descriptive stats with success/failure info from this move
Stats.nTrial = Stats.nTrial + 1;
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

TS = newPsi.TS;
deletedIDs = setdiff( 1:max(keepIDs), keepIDs );
renamedKeepIDs = keepIDs;
for dd = deletedIDs
    renamedKeepIDs( keepIDs > dd ) = renamedKeepIDs( keepIDs > dd ) - 1;
end
for ii = 1:length( TS.obj )
    keepTS = [];
    keepFeatIDs = [];
    for kk = 1:length(keepIDs)
        kID = find( TS.obj(ii).availFeatIDs == keepIDs(kk) );
        if ~isempty( kID )
            keepTS(end+1) = kID;
            keepFeatIDs(end+1) = renamedKeepIDs(kk);
        end
    end
    assert( all( keepFeatIDs == find( F(ii,:) ) ), 'BAD BAD BAD' );
    TS.obj(ii).pi_z = TS.obj(ii).pi_z( keepTS, keepTS );
    TS.obj(ii).availFeatIDs = keepFeatIDs;
    TS.obj(ii).pi_init = ones( 1, length(keepFeatIDs) );
end
theta = newPsi.theta;
switch model.obsModel.type
    case 'Multinomial'
        theta.logp = theta.logp( keepIDs, : );
    case 'Gaussian'
        theta.Mu = theta.Mu( keepIDs, : );
        theta.invSigma = theta.invSigma(:, :, keepIDs );
    case 'AR-Gaussian'
        theta.A = theta.A( :, :, keepIDs );
        theta.invSigma = theta.invSigma(:, :, keepIDs );
    otherwise
        error( 'NEED TO UPDATE THIS!' );
end

end % main function
