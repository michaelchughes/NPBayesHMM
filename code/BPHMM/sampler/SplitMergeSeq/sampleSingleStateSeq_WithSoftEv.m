function [z, logPrQ_z] = sampleSingleStateSeq_WithSoftEv( ii, seqTransM, seqSoftEv, TargetPsi)
% Sample candidate state sequence for a single time series (the ii-th one)
%  given specified transition parameters and soft evidence parameters
%INPUT:
%   ii : integer id of the sequence in question (for using in TargetPsi)
%   seqTransM : transition model struct for seq ii
%           *only* should include features that can appear in state seq.
%            so field eta should be Kii by Kii
%   seqSoftEv : soft evidence matrix for seq ii
%            with dimensions Kii x T
%   TargetPsi : (optional) target configuration to be sampled
%           when this option provided, we *do not* actually sample z
%           but instead compute probability of sampling
%           TargetPsi.stateSeq(ii).z under given Trans/Emission params
%   NOTE: TargetPsi has two crucial fields:
%          .activeFeatIDs specifies the feature IDs in its internal repr.
%          .externalFeatIDs specify the corresponding IDs in seqTransM's
%          representation
%OUTPUT
% z : integer state sequence for time series ii
% logQ : log prob. of sampling sequence z, given input params

logPrQ_z = 0;
T = size( seqSoftEv,2 );

availFeatIDs = seqTransM.availFeatIDs;
if length( availFeatIDs ) == 1
    z = availFeatIDs(1) * ones(1,T);
    return;
end
pi_z = seqTransM.eta;
pi_z = bsxfun( @rdivide, pi_z, sum(pi_z,2)  );

pi_init = 1/length(pi_z)*ones( 1, length(pi_z) );


% Safely convert provided *log* soft evidence matrix
%  into proper probabilities 
%  (up to a prop. constant, so stays in range of CPU)
normC = max( seqSoftEv, [], 1);
seqSoftEv = bsxfun( @minus, seqSoftEv, normC );
likelihood = exp( seqSoftEv );

if ~exist( 'TargetPsi', 'var') || isempty( TargetPsi )
    % ------------------------------------------- Actually sample z(t)
    [js, logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, -1, randi([1 100000]) );
    z = availFeatIDs( js );
else % ----------------------------- Calc prob of moving to Target's z    
    jseq = zeros(1,T);
    for jj = 1:length( availFeatIDs )
        jseq( TargetPsi.stateSeq(ii).z == availFeatIDs(jj) ) = jj;
    end
    [js,logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, jseq, randi([1 100000]) );
    z = availFeatIDs( js );
end

if nargout > 1
    if exist( 'logQ', 'var' )
        logPrQ_z = logPrQ_z + logQ;
    else
        logPrQ_z = logPrQ_z + sum( log( qs ) );
    end
end
assert( ~isnan( logPrQ_z ), 'Bad Q Calculation!' );


end % main function