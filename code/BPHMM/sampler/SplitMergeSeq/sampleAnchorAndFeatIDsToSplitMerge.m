function [anchorIDs, featIDs, PrQfeatIDs] = sampleAnchorAndFeatIDsToSplitMerge( Psi, data, algParams );
% Sample anchor ids and feature ids to perform split/merge move.
% Every split-merge move is defined by:
%   -- two *distinct* sequences (ii,jj), which we call "anchors"
%   -- two features, ki and kj, where F(ii,ki)=1 and F(jj,kj) = 1
% If ki == kj, we propose splitting these into two features,
%   otherwise, we merge ki and kj into a single feature.

F = Psi.F==1;
stateSeq = Psi.stateSeq;
ThetaM = Psi.ThetaM;

% ---------------------------------------------  select anchor sequences
if ~isfield( Psi, 'anchorIDs' ) 
    anchorIDs = randsample( data.N, 2 );
    else
    anchorIDs = Psi.anchorIDs;
end
ii = anchorIDs(1);
jj = anchorIDs(2);

if isfield( Psi, 'activeFeatIDs' )
    ki = Psi.activeFeatIDs(1);
    if length( Psi.activeFeatIDs ) > 1
        kj = Psi.activeFeatIDs(2);
    else
        kj = ki;
    end
end

% ---------------------------------------------  select feature IDs
switch algParams.SM.featSelectDistr
    case 'random'
        qs_ki = F(ii,:);
        if ~exist('ki','var')
            ki = multinomial_single_draw( qs_ki );
        end
        qs_kj = F(jj,:);
        if ~exist('kj','var')
            kj = multinomial_single_draw( qs_kj );
        end
    case 'splitBias'
        qs_ki = F(ii,:);
        if ~exist('ki','var')
            ki = multinomial_single_draw( qs_ki );
        end
        
        delta_ki = false( size(F,2) );
        delta_ki( ki ) = 1;
        
        qs_kj = F(jj,:) .* ~delta_ki;
        qs_kj( ki ) = F(jj,ki)*2*sum( qs_kj );
        if ~exist('kj','var')
            kj = multinomial_single_draw( qs_kj );
        end
    case 'splitBias+margLik'    
        % Build cond distr. kj | ki, jj
        %  based on margLik ratio between ki and kj
        
        qs_ki = F(ii,:);
        if ~exist('ki','var')
            ki = multinomial_single_draw( qs_ki );
        end
        log_qs_kj = -inf( 1, size(F,2)  );
        for kk = find( F(jj,:) )
            if kk == ki
               continue;
            end
            log_qs_kj(kk) = ThetaM.calcMargLikRatio_MergeFeats( data, stateSeq, ki, kk );
        end
        M = max( log_qs_kj );
        if all( isinf(log_qs_kj) )
            qs_kj = zeros(1, size(F,2) );
            qs_kj(ki) = F(jj,ki);
        else
        qs_kj = exp( log_qs_kj - M );
        qs_kj( ki ) = F(jj,ki)*2*sum( qs_kj ); 
        qs_kj = qs_kj ./ sum( qs_kj );
        end
        % Final smoothing: take convex combo of 99% our qs and 1% us
        %   this avoids terrible reverse probabilities
        us = false( 1, size(F,2) );
        us( F(jj,:) ) = 1;
        us = us./sum(us);
        qs_kj = .99*qs_kj + 0.01*us; 
        
        if ~exist('kj','var')
            kj = multinomial_single_draw( qs_kj );    
        end  
end
qs_ki = qs_ki ./ sum( qs_ki );
qs_kj = qs_kj ./ sum( qs_kj );
featIDs = [ki kj];

assert( ~any( isnan(qs_kj) ), 'ERROR: bad numerical calc of feat select distr.');

PrQfeatIDs = qs_ki( ki ) * qs_kj( kj );
