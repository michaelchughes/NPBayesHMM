function [splitPsi, logQ] = getSplitConfig_SeqAllocation( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs, TargetPsi )

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

if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
    TargetPsi.externalFeatIDs = propFeatIDs;
end

% Keep *only* the params for the feature kA, kB being split
%   this is all we update
AXstats  = getXSuffStats( propF, propStateSeq, data_struct, model, [ii jj], [propFeatIDs] );

% Initialize some key sufficient stats
featCounts = sum( propF( anchorObjIDs, propFeatIDs ), 1 );
nObj = length( anchorObjIDs );
counter=0;

for aa = [ otherObjIDs( randperm(length(otherObjIDs)) ) ii jj]
    counter=counter+1;
    
    if aa == ii || aa == jj
        % Make sure to update F suff. stats for anchors
        %  since we've already "seen" them before, unlike other sequences
        featCounts = featCounts - propF(aa, propFeatIDs);
        nObj = nObj - 1;
        
        Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
        AXstats = decXStats( AXstats, Xdelta, model.obsModel );
    end
    
    % Calc Soft Evidence for the current sequence, including both new feats
    curFeatIDs = find( propF(aa,:));
    allPropFeatIDs = union( curFeatIDs, propFeatIDs );
    LL = calcLogCondPrObsGivenTheta( data_struct(aa), myTheta, obsModel, allPropFeatIDs, size(propF,2) );
    
    % Propose new feature assignments for the cur. seq. 
    myFPr = featCounts ./ ( nObj + hyperparams.c );
    [newFaa, logQ_Faa] = sampleFfromRestrictedSet( aa, anchorObjIDs, LL, propFeatIDs, propTS.obj( aa ), myFPr, TargetPsi );
    propF(aa, propFeatIDs) = newFaa;
    
    % Keep track of probability of making this specific proposal
    logQ.F = logQ.F + logQ_Faa;    
        
    % In case we chose not to include one of the proposed new split feats,
    %   we need to remove that feature's transition weights from propTS
    %   and also remove that feature's entries in soft ev. matrix LL
    availFeatIDs = propTS.obj(aa).availFeatIDs;
    kiloc = find( availFeatIDs == propFeatIDs(1) );
    kjloc = find( availFeatIDs == propFeatIDs(2) );
    switch ( newFaa * [10;1] )
        case 01
            keepIDs = [1:kiloc-1 kiloc+1:length(availFeatIDs)];
        case 10
            keepIDs = [1:kjloc-1 kjloc+1:length(availFeatIDs)];
        case 11
            keepIDs = 1:length(availFeatIDs);
    end
    curEta.availFeatIDs = propTS.obj(aa).availFeatIDs( keepIDs );
    curEta.pi_z = propTS.obj(aa).pi_z( keepIDs, keepIDs );    
    activeFeatIDs = union( curFeatIDs, propFeatIDs(newFaa==1) );
    activeLL = LL( activeFeatIDs, : );
    
    
    % Propose new State Seq for current sequence
    %[propStateSeq, logQaa] = sampleZ_RGS(data_struct, propTS, propF, myTheta, obsModel, propStateSeq, aa, TargetPsi);
    [propStateSeq(aa).z, logPrQ_z] = sampleStateSeq_PrecompSoftEv( aa, curEta, activeLL, TargetPsi);
    logQ.z = logQ.z + logPrQ_z;
    
    % Update suff. stats for F sampling
    featCounts = featCounts + propF(aa, propFeatIDs);
    nObj = nObj + 1;
    
    % Update Emission Params for new proposed features
    %   using the newly allocated data from the current sequence
    Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
    AXstats = incXStats( AXstats, Xdelta, model.obsModel );
    myTheta = setThetaToPosteriorMean( myTheta, AXstats, obsModel, [propFeatIDs] );
    
end

% Record final launch state
splitPsi.F = propF;
splitPsi.stateSeq = propStateSeq;
splitPsi.activeFeatIDs = propFeatIDs;
splitPsi.anchorObjIDs  = anchorObjIDs;
end

