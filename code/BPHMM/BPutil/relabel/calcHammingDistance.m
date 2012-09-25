function [HamDist] = calcHammingDistance( TrueSeq, EstSeq )
% Given two sequences, align the EstSeq to TrueSeq,
%   and compute the hamming distance.
% Note that Hamming distance is SYMMETRIC
%   so compare TrueSeq to EstSeq will yield same dist as EstSeq to TrueSeq
% The prescribed order of arguments here only matters
%   in case we have GroundTruth (storing labels in "true_labels" field)
% INPUT
% TrueSeq and EstSeq can be either struct arrays or row vectors
%   struct array: TrueSeq(ii).z is row vector of labels for ii-th seq.          
%   row vector  : 
% EXAMPLE USAGE
% >> SeqA( 1 ).z = [ 1*ones(1,50) 2*ones(1,50) ];
% >> SeqB( 1 ).z = [ 5*ones(1,33) 10*ones(1,33) 11*ones(1,34) ];
% >> calcHammingDistance( SeqA, SeqB )
% Expected Output:  0.33

% ======================================================== Sanitize Input
if isstruct( TrueSeq ) && isstruct( EstSeq )
    % Squash structs into concatenated long sequences
    zTrue = [];
    zEst  = [];    
    for aa = 1:length( TrueSeq )
        if isfield( TrueSeq(aa), 'true_labels' )
            zTrue = [zTrue TrueSeq(aa).true_labels];
        else
            zTrue = [zTrue TrueSeq(aa).z ];
        end
        zEst  = [zEst  EstSeq(aa).z];        
        assert( length(zTrue)==length(zEst), 'ERROR: sequences must be same length' );
    end
else
    zTrue = TrueSeq;
    zEst = EstSeq;
    assert( length(zTrue)==length(zEst), 'ERROR: sequences must be same length' );
end
% Need to make sure the label sets used are compact
%   e.g. instead of zEst = [4 4 5 5 10 10 10]
%    should be totally renamed to be [1 1 2 2 3 3 3]
zTrue0 = map2smallestIntegers( zTrue, max(zTrue));
zEst0  = map2smallestIntegers( zEst,  max(zEst));

% ======================================================== Align Labels
% Construct distance matrix for alignment
[DM uTrue uEst] = buildDistanceMatrix( zTrue0, zEst0);

% Run Munkres (Hungarian) algorithm to get optimal label alignment
%  assignment is a vector with entry for each unique label in zEst0
%  assignment(5) = 3  says that all entries in zEst0 with label 5
%         are best aligned to true label 3
[assignment, cost] = assignmentoptimal( DM );

% Relabel zEst0 to align with zTrue0
Lnew = max(uTrue)+1;
zEstR = zeros( size( zEst0 ) );
for ii=1:length(assignment)
    if assignment(ii)==0
        zEstR( zEst0 == uEst(ii) )= Lnew;
        assignment(ii) = Lnew;
        Lnew = Lnew + 1;
    else
        zEstR( zEst0 == uEst(ii) )= assignment(ii);
    end
end

% ======================================================== Calc Ham Dist
HamDist = sum( zEstR ~= zTrue0 )/length(zTrue0);