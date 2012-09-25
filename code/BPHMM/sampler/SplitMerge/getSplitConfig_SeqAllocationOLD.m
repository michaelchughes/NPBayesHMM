function [splitPsi, logQ] = SMS_LaunchSplit( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs, TargetPsi )

if ~exist( 'TargetPsi', 'var' )
    TargetPsi = [];
end


if isfield( Psi, 'theta' );
    Psi = rmfield( Psi, 'theta');
    Psi = rmfield( Psi, 'TS');
end
nObj = size(Psi.F,1);
logQ = struct( 'F', 0, 'z', 0 );

obsModel = model.obsModel;

ii = anchorObjIDs(1);
jj = anchorObjIDs(2);

kold = featIDs(1);

otherFeatIDs = setdiff( 1:size(Psi.F,2), featIDs );

otherObjIDs = setdiff( find( Psi.F( :, kold ) > 0   )', anchorObjIDs );
allObjIDs = [anchorObjIDs otherObjIDs];

kA = size(Psi.F,2) + 1;
kB = size(Psi.F,2) + 2;
propFeatIDs = [kA kB];

propF = Psi.F;
propF(:, [kold] ) = 0;
propF(ii, [kA kB] ) = [1 0];
propF(jj, [kA kB] ) = [0 1];

% First, draw starting value for theta from sequence ii and jj
propStateSeq = Psi.stateSeq;
propStateSeq( ii ).z(  propStateSeq( ii ).z == kold )  = kA;
propStateSeq( jj ).z(  propStateSeq( jj ).z == kold )  = kB;

% Create deterministic trans params Eta
propTS = getEta_PriorMean( propF, [], hyperparams, model, allObjIDs, propFeatIDs );

% Create deterministic Theta for ALL features: others + [kA kB]
[Xstats] = getXSuffStats( propF, propStateSeq, data_struct, model, 1:nObj, [otherFeatIDs propFeatIDs] );
myTheta = setThetaToPosteriorMean( [], Xstats, obsModel, [otherFeatIDs propFeatIDs] );

% switch obsModel.type
%     case 'Gaussian'
%         myTheta.Mu( propFeatIDs,:)
%     case 'AR-Gaussian';
%                 anchorObjIDs
%
%         Xstats.nObs( end-1:end )
%         myTheta.A( 1:2, 1:2, propFeatIDs )
%         myTheta.invSigma( 1:3, 1:3, propFeatIDs )
% end

if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
    TargetPsi.externalFeatIDs = propFeatIDs;
end
% TO DO: consider an initial estimate for z_ii and z_jj

% Keep *only* the params for the feature kA, kB being split
%   this is all we update
AXstats  = getXSuffStats( propF, propStateSeq, data_struct, model, [ii jj], [propFeatIDs] );

featCounts = sum( propF( anchorObjIDs, propFeatIDs ), 1 );
nObj = length( anchorObjIDs );
counter=0;
%figure; hold on;
for aa = [ otherObjIDs( randperm(length(otherObjIDs)) ) ii jj]
    counter=counter+1;
    if aa == ii || aa == jj
        featCounts = featCounts - propF(aa, propFeatIDs);
        nObj = nObj - 1;
        
        Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
        AXstats = decXStats( AXstats, Xdelta, model.obsModel );
    end
    
    switch aa
        case ii
            restrictedFeatureSet = [1 0; 1 1];
        case jj
            restrictedFeatureSet = [0 1; 1 1];
        otherwise
            restrictedFeatureSet = [0 1; 1 0; 1 1];
    end
    
    curFeatIDs = find( propF(aa,:));
    allPropFeatIDs = union( curFeatIDs, propFeatIDs );
    %logPrObsUnderAllFeatures = compute_likelihood_unnorm( data_struct(aa), myTheta, obsModel, allPropFeatIDs, size(propF,2), 0, aa);
    logPrObsUnderAllFeatures = calcLogCondPrObsGivenTheta( data_struct(aa), myTheta, obsModel, allPropFeatIDs, size(propF,2) );
    logPr = -Inf(1,3);
    for featureVals = restrictedFeatureSet'
        featureVals = featureVals'; % make it row vector again
        featureCombo = (featureVals * [10;1] );
        switch featureCombo
            case 01
                % propPi = TS.obj(ii) with featIDs(1) deleted
                propPi = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(1), 0);
                propPi = propPi.obj(aa);
            case 10
                % propPi = TS.obj(ii) with featIDs(2) deleted
                propPi = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(2), 0);
                propPi = propPi.obj(aa);
            case 11
                propPi = propTS.obj(aa);
        end
        
        %logPrY = calcLogPrAllObsUnderGivenFeatures( data_struct(aa), propPi, logPrObsUnderAllFeatures, obsModel.type );
        logPrY = calcLogMargPrObsSeq(  propPi, logPrObsUnderAllFeatures );
        switch aa
            case anchorObjIDs(1)
                logPrF = featureVals(2) * log( featCounts(2) ) + (1-featureVals(2)).*log( nObj - featCounts(2) +hyperparams.c0 );
            case anchorObjIDs(2)
                logPrF = featureVals(1) * log( featCounts(1) ) + (1-featureVals(1)).*log( nObj - featCounts(1) +hyperparams.c0 );
            otherwise
                logPrF = sum( featureVals .* log( featCounts ) + (1-featureVals).*log( nObj - featCounts +hyperparams.c0 ) );
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
    if isempty( TargetPsi )
        [choice] = multinomial_single_draw( ps );
    else
        targetVals = TargetPsi.F( aa,  TargetPsi.activeFeatIDs );
        %targetValStr = [num2str(targetVals(1)) num2str(targetVals(2))];
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
            propF(aa,propFeatIDs) = [0 1];
            pTS = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(1), 0);
        case 2
            propF(aa,propFeatIDs) = [1 0];
            pTS = buildTransStructForProposal( propTS, hyperparams, model, aa, propFeatIDs(2), 0);
        case 3
            propF(aa,propFeatIDs) = [1 1];
            pTS = propTS;
    end
    
    logQ.F = logQ.F + log( ps(choice)/sum(ps)  );
    
    if isinf( logQ.F )
        assert(1==1,'yes');
    end
    
    % Update suff. stats for F sampling
    featCounts = featCounts + propF(aa, propFeatIDs);
    nObj = nObj + 1;
    
    % Update State Seq
    [propStateSeq, logQaa] = sampleZ_RGS(data_struct, pTS, propF, myTheta, obsModel, propStateSeq, aa, TargetPsi);
    logQ.z = logQ.z + logQaa;
    
    % Update Emission Params
    %curObjStats = getStateSeqSuffStats_C1( propStateSeq, propF, data_struct, model, [aa], propFeatIDs );
    %Ustats.Nkv = Ustats.Nkv + curObjStats.Nkv;
    %myTheta = getTheta_PosteriorMean( myTheta, Ustats, obsModel, propFeatIDs );
    
    
    Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
    AXstats = incXStats( AXstats, Xdelta, model.obsModel );
    myTheta = setThetaToPosteriorMean( myTheta, AXstats, obsModel, [propFeatIDs] );
    
    %     imagesc( Ustats.Nkv );
    %     title( sprintf( 'after %d objects', nObj )  );
    %     pause;
%     if counter==1 || mod( counter, 25 ) == 0
%         switch obsModel.type
%             case 'Gaussian'
%                 myTheta.Mu( propFeatIDs,:)
%             case 'AR-Gaussian';
%                 anchorObjIDs
%                 
%                 Xstats.nObs( end-1:end )
%                 myTheta.A( 1:2, 1:2, propFeatIDs )
%                 %myTheta.invSigma( 1:3, 1:3, propFeatIDs )
%         end
%     end
%     plot( counter, max(logQ.F, -1000), 'ko', 'MarkerSize',12 );
%     plot( counter, logQ.z, 'r.' );
end

% Record final launch state
splitPsi.F = propF;
splitPsi.stateSeq = propStateSeq;
splitPsi.activeFeatIDs = propFeatIDs;
splitPsi.anchorObjIDs  = anchorObjIDs;
end

