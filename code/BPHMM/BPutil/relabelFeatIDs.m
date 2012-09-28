function [Psi, keepIDs] = relabelFeatIDs( Psi, keepIDs )
% Obtain Left-Ordered Form ordering on the features
%   and remove all features not allocated to any time series obj.
% As a result, F is reindexed by swapping columns
%    so it is in Left-Ordered Form and any columns of all zeros are dropped
% Relabelling must also apply to HMM transition and emission params
%    and also to the HMM hidden state sequence and some cached suff. stats
% These are just *relabelled*... their value (in log prob.) does not change

origF = Psi.F;
Psi.F = origF(:, keepIDs);
if isfield( Psi,'cache') && isfield( Psi.cache, 'logSoftEv')
    for ii = 1:size( Psi.F, 1 )
        [Kmax,Tii] = size( Psi.cache.logSoftEv{ii} );
        if max( keepIDs ) <= Kmax
            Psi.cache.logSoftEv{ii} = Psi.cache.logSoftEv{ii}( keepIDs, : );
        else
            kIDs = keepIDs( keepIDs <= Kmax );
            SoftEvTmp = -inf( length(keepIDs), Tii );
            SoftEvTmp( keepIDs <= Kmax, : ) = Psi.cache.logSoftEv{ii}( kIDs, : );
            Psi.cache.logSoftEv{ii} = SoftEvTmp;
        end
    end
end
if isfield( Psi, 'TransM' )
Psi.TransM = Psi.TransM.reallocateFeatIDs( origF, keepIDs );
end
if isfield( Psi, 'ThetaM' )
Psi.ThetaM = Psi.ThetaM.reallocateFeatIDs( keepIDs );
end

if isfield( Psi, 'stateSeq' )
    stateSeq = Psi.stateSeq;
    for ii = 1:size(Psi.F,1)        
        for jj = 1:length( keepIDs )
            stateSeq(ii).z( Psi.stateSeq(ii).z == keepIDs(jj) ) = jj;
        end
    end
    Psi.stateSeq = stateSeq;
end

if isfield( Psi, 'activeFeatIDs' )
   for kk = 1:length( Psi.activeFeatIDs )
        Psi.activeFeatIDs(kk) = find( keepIDs == Psi.activeFeatIDs(kk)  );
   end
end
