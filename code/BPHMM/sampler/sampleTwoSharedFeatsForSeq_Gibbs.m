function [newF, logQ] = sampleTwoSharedFeatsForSeq_Gibbs( aa, propFeatIDs, softEvIN, seqTransMIN, priorPrFeats, anchorObjIDs, TargetPsi);
% Sample feature assignments for sequence aa for the feature pair
%   indicated by propFeatIDs.

switch aa
    case anchorObjIDs(1)
        restrictedFeatureSet = [1 0; 1 1];
    case anchorObjIDs(2)
        restrictedFeatureSet = [0 1; 1 1];
    otherwise
        restrictedFeatureSet = [0 1; 1 0; 1 1];
end

availFeatIDs = seqTransMIN.availFeatIDs;
kiloc = find( availFeatIDs == propFeatIDs(1) );
kjloc = find( availFeatIDs == propFeatIDs(2) );

logPr = -Inf(1,3);
for featureVals = restrictedFeatureSet'
    featureVals = featureVals'; % make it row vector again
    featureCombo = (featureVals * [10;1] );

    switch featureCombo
        case 01
            keepIDs = [1:kiloc-1 kiloc+1:length(availFeatIDs)];
        case 10
            keepIDs = [1:kjloc-1 kjloc+1:length(availFeatIDs)];
        case 11
            keepIDs = 1:length(availFeatIDs);
    end
    seqTransM.availFeatIDs = seqTransMIN.availFeatIDs( keepIDs );
    seqTransM.eta = seqTransMIN.eta( keepIDs, keepIDs );
    
    logPrY = calcLogMargPrObsSeqFAST(  softEvIN(availFeatIDs(keepIDs), :), seqTransM.eta );
    switch aa
        case anchorObjIDs(1)
            logPrF = featureVals(2) * log( priorPrFeats(2) ) + (1-featureVals(2)).*log( 1-priorPrFeats(2) );
        case anchorObjIDs(2)
            logPrF = featureVals(1) * log( priorPrFeats(1) ) + (1-featureVals(1)).*log( 1-priorPrFeats(1) );
        otherwise
            logPrF = sum( featureVals .* log( priorPrFeats ) + (1-featureVals).*log( 1-priorPrFeats ) );
    end
    
    switch featureCombo
        case 01
            logPr(1) = logPrY + logPrF;
        case 10
            logPr(2) = logPrY + logPrF;
        case 11
            logPr(3) = logPrY + logPrF;
    end
end

% --------------------------------------------------- Gibbs sample asgnment
ps = exp( logPr - max(logPr)  );
if isempty( TargetPsi ) || ~exist('TargetPsi', 'var')
    [choice] = multinomial_single_draw( ps );
else
    targetVals = TargetPsi.F( aa,  TargetPsi.activeFeatIDs );
    switch ( targetVals * [10;1] )
        case 01
            choice = 1;
        case 10
            choice = 2;
        case 11
            choice = 3;
    end
end
switch choice
    case 1
        newF = [0 1];
    case 2
        newF = [1 0];
    case 3
        newF = [1 1];
end

logQ = log( ps(choice)/sum(ps)  );

assert( ~isnan( logQ ), 'Bad Q Calculation!' );
