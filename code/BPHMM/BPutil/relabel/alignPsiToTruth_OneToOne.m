function [alignedPsi, HDist] = alignPsiToTruth_OneToOne( Psi, data )
% Align model parameters recovered by inference
%  so that labels correspond exactly to a known *ground truth* state seq
% This uses a ONE to ONE mapping, so for each ground truth state,
%  we find the single best match among estimated states
% This mapping should minimize Hamming distance between the two configs.
%  Note that ThetaM and TransM are *dropped* by this method.
%USAGE
%  [Psi] = alignPsiToTruth_OneToOne( Psi, data )


Zest = horzcat( Psi.stateSeq(:).z );
Ztrue = data.zTrueAll;

unique_true = unique( Ztrue );
unique_est = unique( Zest );
[ZestR, HamDist, AlignedIDs] = mapEstLabels2Truth( Ztrue, Zest);
% AlignedIDs is a 1 x nUniqueEst vector
%   where entry kk says which true label the kk-th feature aligns to

TsCDF = data.aggTs;
F = Psi.F;
alignedSeq = Psi.stateSeq;
for ii = 1:length( alignedSeq );
    alignedSeq(ii).z = ZestR( TsCDF(ii)+1:TsCDF(ii+1)  );
    
    est_z = ZestR( TsCDF(ii)+1:TsCDF(ii+1)  );
    true_z = Ztrue( TsCDF(ii)+1:TsCDF(ii+1) );
    HDist.Ts(ii) = length( true_z );
    HDist.obj(ii) = sum( ( est_z ~= true_z ) )/length( true_z );
end
HDist.total = HamDist;
M = max( size(F,2), length(unique_true) );


if length( unique_est ) <= length( unique_true )
    alignedF = zeros( data.N, M );
    
    alignedF( :, AlignedIDs ) = F(:, unique_est );
    %alignedFfrac = zeros( length(data_struct), M );
    %alignedFfrac( :, AlignedIDs ) = Ffrac(:, unique_est );
    
    % Finally, we need to fill in "bogus/deadbeat" features at the end
    deadIDs = setdiff( 1:size( F,2 ), unique_est );
    if ~isempty( deadIDs )
       alignedF( :, length(unique_est) + (1:length(deadIDs)) ) = F( :, deadIDs ); 
       %alignedFfrac( :, length(unique_est) + (1:length(deadIDs)) ) = Ffrac( :, deadIDs ); 
    end
    
elseif size(F,2) > length( unique_true )
    
    % sortFeatIDs is a 1 x M vector
    %   that when applied as indices to a (expanded) version of F
    %   sorts the cols so that it aligns with Ftrue
    sortFeatIDs = zeros( 1,  M );
    for aa = 1:length( unique_est )
        sortFeatIDs( AlignedIDs(aa) ) = unique_est(aa);
    end
    
    % If we have too many estimated features ( more than # true states)
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
    %alignedFfrac = Ffrac( :, sortFeatIDs );
end

alignedPsi.F = alignedF;
alignedPsi.stateSeq = alignedSeq;

%alignedPsi.ThetaM = Psi.ThetaM.reallocateFeatIDs(  sortFeatIDs );
%alignedPsi.TransM = Psi.TransM.reallocateFeatIDs(  sortFeatIDs );


end

