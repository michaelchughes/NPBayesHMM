function [relabeled_est_labels, hamming_dist, assignment] = mapEstLabels2Truth(true_labels,est_labels)
% Function uses Munkres algorithm to assign best correspondence in labels
%   between a "ground truth" label sequence and an estimated sequence

small_labels = map2smallestIntegers(true_labels,max(true_labels));
assert( all( unique(small_labels) == unique(true_labels) ), 'BAD TRUE LABELS' );

[DM unique_true unique_est] = buildDistanceMatrix(true_labels,est_labels);

if isempty(setdiff(0,unique_true))
    error('need to shift true_labels')
elseif isempty(setdiff(0,unique_est))
    error('need to shift est_labels')
end

[assignment, cost] = assignmentoptimal(DM);

new_label_count=max(unique_true)+1;

relabeled_est_labels = est_labels;
for ii=1:length(assignment)
    if assignment(ii)==0
        relabeled_est_labels( est_labels==unique_est(ii) )=new_label_count;
        assignment(ii) = new_label_count;
        new_label_count = new_label_count + 1;
    else
        relabeled_est_labels( est_labels==unique_est(ii) )=assignment(ii);
    end
end

hamming_dist = sum(relabeled_est_labels~=true_labels)/length(true_labels);