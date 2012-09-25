function [propF, propTS, logQ_RGS] = sampleFeatures_RGS( F, TS, theta, data_struct, hyperparams, model, anchorObjIDs, propFeatIDs, TargetPsi )
% Sample binary feature asgn matrix F  under RESTRICTED settings
%   which means we
%       (1) only update features for sequences in list "objIDs"
%       (2) only consider asgns where at least one of "featIDs" is positive
% Executed under two circumstances:
%   (I)  actually sample from posterior on F over [0 1], [1 0], [1 1]
%   (II) do not sample at all, but compute prob.
%            for moving from current F to TargetPsi.F with fixed TS,theta

propF = F;
logQ_RGS = 0;
obsModel = model.obsModel;

activeObjIDs = find( sum(F(:, propFeatIDs), 2) > 0  )';

% Expand TS to make sure all objects possess both feature ids
propTS = expandTStoIncludeFeatIDs( F, TS, hyperparams, activeObjIDs, propFeatIDs );

% Initialize some key sufficient stats
featCounts = sum( F( activeObjIDs, propFeatIDs ), 1 );
nObj = size(F,1);

didWarn=0;
for aa = activeObjIDs

     % Make sure to update F suff. stats
     featCounts = featCounts - propF(aa, propFeatIDs);
     nObj = nObj - 1;
    
    % Calc Soft Evidence for the current sequence, including both new feats
    curFeatIDs = find( propF(aa,:));
    allPropFeatIDs = union( curFeatIDs, propFeatIDs );
    LL = calcLogCondPrObsGivenTheta( data_struct(aa), theta, obsModel, allPropFeatIDs, size(propF,2) );
    
    % Propose new feature assignments for the cur. seq. 
    myFPr = featCounts ./ ( nObj + hyperparams.c );
    [newFaa, logQ_Faa] = sampleFfromRestrictedSet( aa, anchorObjIDs, LL, propFeatIDs, propTS.obj( aa ), myFPr, TargetPsi );
    propF(aa, propFeatIDs) = newFaa;
    
    % Keep track of probability of making this specific proposal
    logQ_RGS = logQ_RGS + logQ_Faa;    
    if ~isinf( logQ_RGS ) && ~didWarn
        %fprintf( 'WARNING: Bad Calculation of Proposal Prob. for F\n' );
        didWarn=1;
    end
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
    propTS.obj(aa).availFeatIDs = propTS.obj(aa).availFeatIDs( keepIDs );
    propTS.obj(aa).pi_z = propTS.obj(aa).pi_z( keepIDs, keepIDs );  
    propTS.obj(aa).pi_init = propTS.obj(aa).pi_init( keepIDs );

    % Update suff. stats for F sampling
    featCounts = featCounts + propF(aa, propFeatIDs);
    nObj = nObj + 1;

end