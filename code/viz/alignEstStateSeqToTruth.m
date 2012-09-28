function [zEstR, TrueAlphaNames, EstAlphaNames] = alignEstStateSeqToTruth( zEst, zTrue, alignP )
% Determines best alignment for each estimated state to ground truth
% Each true state is aligned to potentially several estimated states
%   with the option that no est state aligns properly
% Each estimated state is aligned to at most 1 true state
%   with option that if insufficient overlap exists, we mark it as "unmatched"
% This alignment is good for *plotting* the stateSeq
%   but not good at all for computing Hamming distance

% To qualify as aligned, an estimated state must have at least
%   OVERLAP_THR fraction of its labels matched to the ground truth sequence
alignP.OVERLAP_THR = 0.5;
alignP.nMISMATCH   = 4;
alignP.NO_MATCH_LABEL = 0;

NO_MATCH_LABEL = alignP.NO_MATCH_LABEL;

% Compute the best true label match of each estimated label
[A, uTrue, uEst] = buildAlignmentMatrix( zTrue, zEst );
zEstR = zeros( size( zEst ) );
uEstR = zeros( size( uEst ) );

% Identify best true state for each est. state
%   bestMatches(kk) gives ID of true state best aligned to kk-th state
[~, bestMatches] = max( A, [], 2 );
MatchCounts = histc( bestMatches, uTrue );

JITTER_AMT = round( 100/(1+max(MatchCounts) ) ) / 100;
jitter = zeros( size( uTrue ) );

% Loop over all est states,
%  and assign to a specific true state, or  "unmatched"
for ue = 1:length( uEst )
    ut = bestMatches(ue);
    if A( ue, ut ) > alignP.OVERLAP_THR
        if sum( zEstR == uTrue( ut ) ) > 0
            % if we've already found a match for this true state
            %   we need to add some jitter
            jitter(ut) = jitter(ut) - JITTER_AMT;
        end
        zEstR( zEst == uEst(ue) ) = uTrue(ut) +jitter(ut);
        uEstR( ue ) = uTrue(ut) +jitter(ut);
    else
        % We need to mark this state as "unmatched" to any ground truth
        NO_MATCH_LABEL = NO_MATCH_LABEL - 1;
        zEstR( zEst == uEst(ue) ) = NO_MATCH_LABEL;
        uEstR( ue ) = NO_MATCH_LABEL;
    end
end

if length( uTrue ) > 26
    fprintf( 'WARNING: probably symbol clashing about to happen... Need more letters!\n');
end
for aa = 1:length( uTrue )
    % TO DO: worry about going over 26 letters for the alphabet!
    TrueAlphaNames{aa} = char(65 + aa-1);
end

JitterCounts = zeros( size( MatchCounts) );
%for ue = 1:length( uEstR )
    %ueOrig = find( uEst == unique( zEst( zEstR == uEstR(ue) ) ) );
    %ut = bestMatches( ueOrig );
for ue = 1:length( uEst )
    ut = bestMatches( ue );
    if uEstR(ue) < 0
        EstAlphaNames{ue} = char(123 + uEstR(ue ));  % select z,y,x,w,...
    elseif MatchCounts(ut) > 1
        JitterCounts(ut) = JitterCounts(ut)+1;
        jC = JitterCounts(ut);
        EstAlphaNames{ue} = [TrueAlphaNames{ut} '_' num2str(jC)];
    else
        EstAlphaNames{ue} = TrueAlphaNames{ut};
    end
end