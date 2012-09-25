function [re_true_labels, re_est_labels] = map2ClosestGroundTruth( true_labels, est_labels )
% Reassign values to est_labels and true_labels
%   in order to make their correspondences more easily seen

true_labels = map2smallestIntegers(true_labels,max(true_labels));
%est_labels = map2smallestIntegers(est_labels,max(est_labels));

% DM(jj,kk) gives distance from est label jj to true_label kk
%   smaller value = better match
[DM unique_true unique_est] = buildDistanceMatrix(true_labels,est_labels);

if isempty(setdiff(0,unique_true))
    error('need to shift true_labels')
elseif isempty(setdiff(0,unique_est))
    error('need to shift est_labels')
end

[~, best_match] = min( DM, [], 2 );

zz = 1;
zStep = 1;
zSKIP = 5;
for kk = 1:length( unique_true )
    re_true_labels( true_labels == unique_true(kk) ) = zz;
    
    nMatch = 0;
    for jj = 1:length( unique_est )
       if best_match(jj) ~= kk
           continue;
       end
       nMatch = nMatch + 1;
       re_est_labels( est_labels == unique_est(jj) ) = zz + (nMatch-1)*zStep;
    end
    
    %est_matches = find( best_match == kk );
    %for ii = 1:length( est_matches )
    %    jj = est_matches(ii);
    %    re_est_labels( est_labels == jj ) = zz + (ii-1)*zStep;
    %end
    
    zz = zz + (nMatch+zSKIP)*zStep;
end

badIDs = find( re_est_labels == 0 );
assert( isempty(badIDs), 'did not label properly' );

%zMax = max( [re_est_labels re_true_labels] );
%BestColorMap = jet( zMax );