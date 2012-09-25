function [alignedF, alignedFfrac, alignedSeq, HDist, nTrue, sortFeatIDs] = mapStateSeqLabelsToGroundTruth( F, Ffrac, stateSeq, data_struct )
% Align feature matrix F and state sequence recovered by inference
%  so that labels correspond to a known *ground truth* state sequence
% Output aligned F and state sequence
%    both correspond to the same ground truth sequences
% Also returns the Hamming distance between the aligned seq and truth
%    in struct format
% INPUT:
%  F     : N x K matrix
%           binary feature matrix
%  Ffrac : N x K matrix
%           real matrix, indicates average over number of posterior samps
%           for F
%  stateSeq : Nx1 struct
%           stateSeq(ii).z gives state labels for sequence ii
%  data_struct : Nx1 struct
%           data_struct(ii).true_labels gives TRUE labels for sequence ii
for ii = 1:length( stateSeq )
    Ts(ii) = length( stateSeq(ii).z );
end
est_seq = zeros( 1, sum( Ts ) );
true_seq = zeros( 1, sum(Ts)  );
TsCDF = [0 cumsum(Ts)];
for ii = 1:length( data_struct )
    Ttrue =  data_struct(ii).T;
    if Ttrue == Ts(ii)
        true_seq(  TsCDF(ii)+1:TsCDF(ii+1) ) = data_struct(ii).true_labels;
    elseif Ttrue == Ts(ii) + 1  % handle auto-regressive case!
        true_seq(  TsCDF(ii)+1:TsCDF(ii+1) ) = data_struct(ii).true_labels(2:end);
    else
        error( 'Mismatch in lengths between true and est sequences!' );
    end
    est_seq(   TsCDF(ii)+1:TsCDF(ii+1) ) = stateSeq(ii).z;
end

unique_true = unique( true_seq );
nTrue = length( unique_true );
[new_est_seq, HamDist, AlignedIDs] = mapEstLabels2Truth(true_seq,est_seq);
% AlignedIDs is a 1 x nUniqueEst vector
%   where entry kk says which true label the kk-th feature aligns to


alignedSeq = stateSeq;
for ii = 1:length( alignedSeq );
    alignedSeq(ii).z = new_est_seq( TsCDF(ii)+1:TsCDF(ii+1)  );
    
    est_z = new_est_seq( TsCDF(ii)+1:TsCDF(ii+1)  );
    true_z = true_seq( TsCDF(ii)+1:TsCDF(ii+1) );
    HDist.Ts(ii) = length( true_z );
    HDist.obj(ii) = sum( ( est_z ~= true_z ) )/length( true_z );
end
HDist.total = HamDist;
M = max( size(F,2), length(unique_true) );
unique_est = unique( est_seq );


if length( unique_est ) <= length( unique_true )
    alignedF = zeros( length(data_struct), M );
    alignedFfrac = zeros( length(data_struct), M );
    
    alignedF( :, AlignedIDs ) = F(:, unique_est );
    alignedFfrac( :, AlignedIDs ) = Ffrac(:, unique_est );
    
    % Finally, we need to fill in "bogus/deadbeat" features at the end
    deadIDs = setdiff( 1:size( F,2 ), unique_est );
    if ~isempty( deadIDs )
       alignedF( :, length(unique_est) + (1:length(deadIDs)) ) = F( :, deadIDs ); 
       alignedFfrac( :, length(unique_est) + (1:length(deadIDs)) ) = Ffrac( :, deadIDs ); 
    end
    
elseif size(F,2) > length( unique_true )
    
    % sortFeatIDs is a 1 x M vector
    %   that when applied as indices to a (expanded) version of F
    %   sorts the cols so that it aligns with Ftrue
    sortFeatIDs = zeros( 1,  M );
    for aa = 1:length( unique_est )
        sortFeatIDs( AlignedIDs(aa) ) = unique_est(aa);
    end
    
    % If we have too many unique features ( more than # true states)
    %   then we need to fill in "zero" entries so that unaligned est states
    %    appear at the end of the aligned matrix (with high feat IDs )
    Ks = setdiff( 1:size(F,2), unique_est );
    if ~isempty( Ks )
        kk=1;
        for aa = 1:size(F,2)
            if sortFeatIDs(aa) == 0
                sortFeatIDs( aa ) = Ks(kk);
                kk = kk + 1;
            end
        end
    end
    
    alignedF = F( :, sortFeatIDs );
    alignedFfrac = Ffrac( :, sortFeatIDs );
end


end



% If we have too many true features (more than # est states )
%    then we need to add in a zero column for every "zero" in sortFeatIDs
%     Ks = setdiff( 1:M, unique_est );
%     emptyIDs = find( sortFeatIDs == 0 );
%     K  = size(F,2);
%     kk = 1;
%     for aa = 1:M
%         if sortFeatIDs(aa) == 0
%             sortFeatIDs( aa ) = Ks(kk );
%             kk = kk+1;
%             %K = K + 1;
%             %sortFeatIDs( aa ) = K;
%         end
%     end
%F = [F zeros(length(data_struct), length(emptyIDs) )];
%Ffrac = [Ffrac zeros( length(data_struct), length(emptyIDs) ) ];
