function [Psi, Stats] = sampleSharedFeats( Psi, data, seqIDs, featIDs )
% Sample all "shared" entries in binary feature matrix F
% Recall that "shared" entries imply sequences ii and feature IDs kk
%   such that the column F(:,kk) has >= one other ON entry other than ii
% This update skips over "unique" entries in the matrix F.
% Fast metropolis-hastings proposals are used, and calculated 
%   log likelihoods cached for efficient access by later sampling methods
%INPUT
%  Psi : input model config... it's field "F" is the matrix to update
%  data : SeqData object, holds observations for all sequences
%  seqIDs : (optional) defines which sequence IDs to sample
%  featIDs : (optional) defines which feature IDs to sample (other ignored)
%OUTPUT
%  Psi : resulting model configuration
%  Stats : records how often each type of proposal was accepted/rejected
%          Stats.C is a 2x2 matrix, counting MH proposal results as follows
%              C=  [ # trans 0 to 0      # trans 0 to 1; 
%                    # trans 1 to 0      # trans 1 to 1]

% -------------------- Unpack
F = Psi.F == 1; % force to type logical
[N,K] = size(F);
if ~exist( 'featIDs', 'var' )
    featIDs = 1:K;
end
if ~exist( 'seqIDs', 'var' )
    seqIDs = 1:N;
end
Stats = struct('C', zeros(2,2) );

% ------------------- Precompute Proposed Eta + Soft Evidence 
PropEta = cell(N,1);
logSoftEv = cell(N,1);
logMargPrObs = -Inf(N,1);
for ii = seqIDs
    PropEta{ii} = Psi.TransM.sampleEtaProposal_Shared( ii );
    logSoftEv{ii} = Psi.ThetaM.calcLogSoftEv( ii, data );
    ks = F(ii,:);
    logMargPrObs(ii) = calcLogMargPrObsSeqFAST( logSoftEv{ii}(ks,:), PropEta{ii}(ks, ks) );
end

% Suff statistics for IBP sampling of F
featCounts = sum( F, 1 );

% ------------------------------------------------- Sample each entry in F
for ii = seqIDs( randperm( length(seqIDs)  ) )
    for kk = featIDs( randperm( length(featIDs) ) )
        
        if featCounts(kk) == 1 && F(ii,kk) == 1
            % Skip updates to unique features!
            continue;
        end
        
        featCounts(kk) = featCounts(kk) - F(ii,kk );
        
        propks = F(ii,:) == 1;
        if F(ii, kk )
           oldF = 1;
           propks(kk) = false;
        else
           oldF = 0;
           propks(kk) = true;
        end
        propEta = PropEta{ii}( propks, propks );
        propLL = logSoftEv{ii}( propks, : );
        
        [ F(ii,kk), logMargPrObs(ii) ] = sampleSingleFeatEntry_SharedMH( F(ii,kk), featCounts(kk), N, Psi.bpM.c, propEta, propLL, logMargPrObs(ii) );

        featCounts(kk) = featCounts(kk) + F(ii,kk );

        Stats.C( oldF+1, F(ii,kk)+1) = Stats.C( oldF+1, F(ii,kk)+1) + 1;
    end
end

% -------------------- Repack
Psi.F = F;
Psi.TransM = Psi.TransM.updateAllEta( F, PropEta, seqIDs );
Psi.cache.logSoftEv = logSoftEv;
Psi.cache.logMargPrObs = logMargPrObs;
Psi = reallocateFeatIDs( Psi );
end % main function
