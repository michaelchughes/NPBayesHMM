function [newF, logQaa] = sampleFfromRestrictedSet( aa, anchorObjIDs, LL, propFeatIDs, propPiIN, myFPr, TargetPsi )
% Sample feature assignments to the feature pair indicated by propFeatIDs

switch aa
    case anchorObjIDs(1)
        restrictedFeatureSet = [1 0; 1 1];
    case anchorObjIDs(2)
        restrictedFeatureSet = [0 1; 1 1];
    otherwise
        restrictedFeatureSet = [0 1; 1 0; 1 1];
end

availFeatIDs = propPiIN.availFeatIDs;
kiloc = find( availFeatIDs == propFeatIDs(1) );
kjloc = find( availFeatIDs == propFeatIDs(2) );

logPr = -Inf(1,3);
for featureVals = restrictedFeatureSet'
    featureVals = featureVals'; % make it row vector again
    featureCombo = (featureVals * [10;1] );

    switch featureCombo
        case 01
            % propPi = TS.obj(ii) with featIDs(1) deleted
            %propPi = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(1), 0);
            %propPi = propPi.obj(aa);
            %propPi = propPiIN;
            keepIDs = [1:kiloc-1 kiloc+1:length(availFeatIDs)];
        case 10
            keepIDs = [1:kjloc-1 kjloc+1:length(availFeatIDs)];
            % propPi = TS.obj(ii) with featIDs(2) deleted
            %propPi = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(2), 0);
            %propPi = propPi.obj(aa);
        case 11
            keepIDs = 1:length(availFeatIDs);
            %propPi = propPiIN; % propTS.obj(aa);
    end
    propPi.availFeatIDs = propPiIN.availFeatIDs( keepIDs );
    propPi.pi_z = propPiIN.pi_z( keepIDs, keepIDs );
    propPi.pi_init = propPiIN.pi_init( keepIDs );
    
    logPrY = calcLogMargPrObsSeqFAST(  LL( propPi.availFeatIDs, :), propPi.pi_z );
    switch aa
        case anchorObjIDs(1)
            logPrF = featureVals(2) * log( myFPr(2) ) + (1-featureVals(2)).*log( 1-myFPr(2) );
        case anchorObjIDs(2)
            logPrF = featureVals(1) * log( myFPr(1) ) + (1-featureVals(1)).*log( 1-myFPr(1) );
        otherwise
            logPrF = sum( featureVals .* log( myFPr ) + (1-featureVals).*log( 1-myFPr ) );
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

logQaa = log( ps(choice)/sum(ps)  );

%assert( ~isinf( logQaa ), 'Bad Calculation of Proposal Prob.' );
%if isinf( logQaa )
%    fprintf( 'WARNING: Forced transition highly improbable (-Inf)\n' );
%end

