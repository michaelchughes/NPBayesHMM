function [z, logPrQ_z] = sampleStateSeq_PrecompSoftEv( ii, Etaii, LL, TargetPsi)
% Sample hidden states "z" for one specific sequence
%INPUT
%  ii : integer ID of which sequence to sample
%  Etaii : struct for ii's eta transition weights (not normalized)
%  LL    : Kii x T   matrix of soft evidence
%  TargetPsi : (optional) desired target values for this state seq
%               are recorded in TargetPsi.stateSeq(ii).z
%             when present, we don't actually sample "z", but calculate
%               the probability of achieving this target configuration

logPrQ_z = 0;
T = size( LL,2 );

if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
    % Translate from TargetPsi's feature IDs to
    %    corresponding IDs in F, theta, prevStateSeq, etc.
    z = TargetPsi.stateSeq(ii).z;
    for aa = 1:length( TargetPsi.activeFeatIDs )
        z( z == TargetPsi.activeFeatIDs(aa) ) = TargetPsi.externalFeatIDs(aa);
    end
    TargetPsi.stateSeq(ii).z = z;
end

availFeatIDs = Etaii.availFeatIDs;
if length( availFeatIDs ) == 1
    z = availFeatIDs(1) * ones(1,T);
    return;
end
pi_z = Etaii.pi_z;
pi_z = bsxfun( @rdivide, pi_z, sum(pi_z,2)  );

pi_init = 1/length(pi_z)*ones( 1, length(pi_z) );


% Safely convert *log* soft evidence matrix into proper probabilities 
%  (up to a prop. constant, so stays in range of CPU)
normC = max( LL, [], 1);
LL = bsxfun( @minus, LL, normC );
likelihood = exp( LL );


MEXSEED = randomseed();
randomseed( MEXSEED+1 );
if ~exist( 'TargetPsi', 'var') || isempty( TargetPsi )
    % ------------------------------------------- Actually sample z(t)
    [js, logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, -1, MEXSEED );
    z = availFeatIDs( js );
else % ----------------------------- Calc prob of moving to Target's z    
    jseq = zeros(1,T);
    for jj = 1:length( availFeatIDs )
        jseq( TargetPsi.stateSeq(ii).z == availFeatIDs(jj) ) = jj;
    end
    [js,logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, jseq, MEXSEED );
    z = availFeatIDs( js );
end

if nargout > 1
    if exist( 'logQ', 'var' )
        logPrQ_z = logPrQ_z + logQ;
    else
        logPrQ_z = logPrQ_z + sum( log( qs ) );
    end
end

end % main function