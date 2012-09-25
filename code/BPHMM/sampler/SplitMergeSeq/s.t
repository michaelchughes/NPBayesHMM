function [anchorIDs, featIDs, PrQfeatIDs] = sampleAnchorAndFeatIDsToSplitMerge( Psi, data, algParams );
% Sample anchor ids and feature ids to perform split/merge move.
% Every split-merge move is defined by:
%   -- two *distinct* sequences (ii,jj), which we call "anchors"
%   -- two features, ki and kj, where F(ii,ki)=1 and F(jj,kj) = 1
% If ki == kj, we propose splitting these into two features,
%   otherwise, we merge ki and kj into a single feature.

F = Psi.F;
stateSeq = Psi.stateSeq;
ThetaM = Psi.ThetaM;

% ---------------------------------------------  select anchor sequences
if ~isfield( Psi, 'anchorIDs' ) 
    anchorIDs = randsample( data.N, 2 );
end
ii = anchorIDs(1);
jj = anchorIDs(2);

if isfield( Psi, 'activeFeatIDs' )
    [ki kj] = Psi.activeFeatIDs;
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
            log_qs_kj(kk) = ThetaM.calcMargLikRatio_MergeFeats( stateSeq, data, [ki kk] );
        end
        M = max( log_qs_kj( log_qs_kj ~= 0 ) );
        qs_kj = exp( log_qs_kj - M );
        qs_kj( ki ) = F(jj,ki)*2*sum( qs_kj );  
                if ~exist('kj','var')
        kj = multinomial_single_draw( qs_kj );    
        end  
end

featIDs = [ki kj];

PrQfeatIDs = qs_ki( ki ) * qs_kj( kj );
