function [mergePsi, logQ] = getMergeConfig_SeqAllocation( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs, TargetPsi )
% Merge Move launch
%   MERGE all blocks assigned to either featIDs into single feature kmerge
if ~exist( 'TargetPsi', 'var' )
    TargetPsi = [];
end

if isfield( Psi, 'theta' );
    Psi = rmfield( Psi, 'theta');
    Psi = rmfield( Psi, 'TS');
end

logQ.F   = 0;
logQ.z   = 0;
obsModel = model.obsModel;

F = Psi.F;
nObj = size(F,1);
ii = anchorObjIDs(1);
jj = anchorObjIDs(2);

allObjIDs = find( sum( F(:,featIDs), 2 ) > 0 )';
otherObjIDs = setdiff( allObjIDs, anchorObjIDs );
otherFeatIDs = setdiff( 1:size(F,2), featIDs );

propF = F;
propF(:, featIDs) = 0;
kmerge = size(propF,2)+1;

propFeatIDs = [kmerge];
propF(allObjIDs, kmerge) = 1;

propStateSeq = Psi.stateSeq;
for aa = anchorObjIDs
    aIDs =  propStateSeq(aa).z == featIDs(1);
    bIDs =  propStateSeq(aa).z == featIDs(2);
    propStateSeq(aa).z( aIDs | bIDs ) = kmerge;
end
propTS = getEta_PriorMean( propF, [], hyperparams, model, allObjIDs, propFeatIDs );

% Ustats  = getStateSeqSuffStats_C1( propStateSeq, propF, data_struct, model, 1:nObj, [otherFeatIDs propFeatIDs] );
% myTheta.logp = zeros( size(Ustats.Nkv)  );
% myTheta = getTheta_PosteriorMean( myTheta, Ustats, obsModel, [otherFeatIDs propFeatIDs] );

[Xstats] = getXSuffStats( propF, propStateSeq, data_struct, model, 1:nObj, [otherFeatIDs propFeatIDs] );
myTheta = setThetaToPosteriorMean( [], Xstats, obsModel, [otherFeatIDs propFeatIDs] );

% Keep *only* the params for the feature kM being merged
%   this is all we update
AXstats  = getXSuffStats( propF, propStateSeq, data_struct, model, [ii jj], [propFeatIDs] );


if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
    TargetPsi.externalFeatIDs = kmerge;
end
for aa = [ otherObjIDs( randperm( length(otherObjIDs) ) ) anchorObjIDs]
    if aa == ii || aa == jj        
        Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
        AXstats = decXStats( AXstats, Xdelta, model.obsModel );    
    end

    % Construct soft evidence
    activeFeatIDs = find( propF(aa,:)==1 );
    LL = calcLogCondPrObsGivenTheta( data_struct(aa), myTheta, obsModel, activeFeatIDs, size(propF,2) );
    LL = LL( activeFeatIDs, :);
    % Propose new State Seq. for current item
    [propStateSeq(aa).z, logPrQ_z] = sampleStateSeq_PrecompSoftEv( aa, propTS.obj(aa), LL, TargetPsi);    
    %[propStateSeq, logQaa] = sampleZ_RGS(data_struct, propTS, propF, myTheta, model.obsModel, propStateSeq, aa, TargetPsi);
    logQ.z = logQ.z + logPrQ_z;
    
    % Update Emission Params for new proposed features
    %   using the newly allocated data from the current sequence
    Xdelta = getXSuffStats( propF, propStateSeq, data_struct, model, aa, propFeatIDs );
    AXstats = incXStats( AXstats, Xdelta, model.obsModel );
    myTheta = setThetaToPosteriorMean( myTheta, AXstats, obsModel, propFeatIDs );
end

mergePsi.F = propF;
mergePsi.stateSeq = propStateSeq;
mergePsi.activeFeatIDs = propFeatIDs;
mergePsi.anchorObjIDs  = anchorObjIDs;