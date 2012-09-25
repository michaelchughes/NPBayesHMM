function ps = getSplitMergeFeatChoiceProbsBetter( F, stateSeq, data_struct, model, jj, ki, Xstats )

availFeatIDs = find( F(jj,:) );
distinctFeatIDs = setdiff( availFeatIDs, ki );

ps = zeros( 1, size(F,2) );

if ~isempty( distinctFeatIDs )
    if ~exist( 'Xstats', 'var' )
        Kstats = getXSuffStats( F, stateSeq, data_struct, model, 1:size(F,1),  ki );
        Dstats = getXSuffStats( F, stateSeq, data_struct, model, 1:size(F,1),  distinctFeatIDs );
    else
        Kstats = selectFromXSuffStats( Xstats, model, ki );
        Dstats = selectFromXSuffStats( Xstats, model, distinctFeatIDs );
    end
    
    sepLogPr = zeros(1, length(distinctFeatIDs) );
    comboLogPr = zeros(1, length(distinctFeatIDs) );
    for rr = 1:length( distinctFeatIDs )
        
        switch model.obsModel.type
            case 'Multinomial'
                SEPARATE = Kstats;
                SEPARATE.Nkv(2,:) = Dstats.Nkv(rr,:);
                sepLogPr(rr) = calcMargLogPrData_MultinomialDirichlet( SEPARATE.Nkv, model.obsModel.params.lambda', 1 );
                
                COMBO = struct();
                COMBO.Nkv = sum( SEPARATE.Nkv, 1);
                comboLogPr(rr) = calcMargLogPrData_MultinomialDirichlet( COMBO.Nkv, model.obsModel.params.lambda', 1 );
                
            case 'Gaussian'
                SEPARATE = Kstats;
                SEPARATE.nObs(2) = Dstats.nObs(rr);
                SEPARATE.Mu(2,:) = Dstats.Mu(rr,:);
                SEPARATE.Sigma(:,:,2) = Dstats.Sigma(:,:,rr);
                sepLogPr(rr) = calcMargLogPrData_GaussianNormalInvWishart( SEPARATE, model.obsModel.params );

                COMBO.nObs = sum( SEPARATE.nObs );
                if COMBO.nObs > 0
                    COMBO.Mu = SEPARATE.nObs(1)*SEPARATE.Mu(1,:) + SEPARATE.nObs(2)*SEPARATE.Mu(2,:);
                    COMBO.Mu = COMBO.Mu / COMBO.nObs;
                    
                    [K D] = size( COMBO.Mu );
                    COMBO.Sigma = zeros( D, D, K );
                    % Can't just add in the "Sigma fields", since they depend
                    % on subtracting the just computed Mu field
                    for ii = 1:length( data_struct )
                        comboINDS = stateSeq(ii).z == distinctFeatIDs(rr) | stateSeq(ii).z == ki;
                        if sum( comboINDS ) > 0
                            obsDiff = bsxfun( @minus, data_struct(ii).obs( comboINDS, : ), COMBO.Mu );
                            COMBO.Sigma = COMBO.Sigma + obsDiff'*obsDiff;
                        end
                    end
                else
                    DD = length( SEPARATE.Mu(1,:) );
                    COMBO.Mu = zeros( 1, DD );
                    COMBO.Sigma = zeros( DD, DD );
                end
                comboLogPr(rr) = calcMargLogPrData_GaussianNormalInvWishart( COMBO, model.obsModel.params );
            case 'AR-Gaussian'
                SEPARATE = Kstats;
                SEPARATE.nObs(2) = Dstats.nObs( rr );
                SEPARATE.XX(:,:,2) = Dstats.XX(:,:,rr);
                SEPARATE.XY(:,:,2) = Dstats.XY(:,:,rr);
                SEPARATE.YY(:,:,2) = Dstats.YY(:,:,rr);
                %             SEPARATE.nObs = Kstats.XX
                %             SEPARATE.XX = Xstats.XX(:, :, featIDs );
                %             SEPARATE.XY = Xstats.XY(:, :, featIDs );
                %             SEPARATE.YY = Xstats.YY(:, :, featIDs );
                sepLogPr(rr) = calcMargLogPrData_ARMatrixNormalInvWishart( SEPARATE, model.obsModel.params );
                
                COMBO.nObs = sum( SEPARATE.nObs );
                COMBO.XX = sum( SEPARATE.XX, 3 );
                COMBO.XY = sum( SEPARATE.XY, 3);
                COMBO.YY = sum( SEPARATE.YY, 3);
                comboLogPr(rr) = calcMargLogPrData_ARMatrixNormalInvWishart( COMBO, model.obsModel.params );
        end
    end

    logps = comboLogPr - sepLogPr;
    assert( all(~isnan(logps)), 'ERROR: bad calc of log prob' );
    M = max( logps( sepLogPr ~= 0 ) );
    ps(distinctFeatIDs) = exp( logps - M );
end

if F( jj, ki ) 
    if ~isempty( distinctFeatIDs )        
        ps( ki ) = 2*( sum( ps( [1:ki-1 ki+1:end] ) ) );
    else
        ps(ki) = 1;
    end
else
    ps(ki)=0;
end

ps = ps./sum(ps);

us = zeros( 1, size(F,2) );
us( availFeatIDs ) = 1;
us = us./sum(us);

% Final smoothing: take convex combo of 99% our ps and 1% us (to avoid terrible reverse probabilities)
ps = .99*ps + 0.01*us;

% EPS = 1e-9;
% switch model.obsModel.type
%     case 'Multinomial'
%         pHat = exp( myTheta.logp )';
%         pHat = bsxfun( @rdivide, pHat, sum(pHat,1) );
%         D = calcChiSqDistanceBetter( pHat(:,ki), pHat );
%     case 'Gaussian'
%         if size( myTheta.Mu, 1) > 1
%             D = pdist( myTheta.Mu, 'euclidean');        
%             D = squareform(D);
%             D = D(ki,:);
%         else
%             D = 0; % Default distance from self is zero.
%         end
%         D = D./(EPS+max(max(D)) ); % Ensure D varies from 0 to 1. 
%         
%     case 'AR-Gaussian'
%         [D DR K] = size( myTheta.A );
%         
%         %Avec = -100*ones(K,D);
%         Avec = zeros(K, D*DR);
%         for kk = relFeatIDs
%             %Sig = myTheta.invSigma(:,:,kk) \ eye(D);
%             %Avec(kk,:) = diag( Sig );
%             Avec(kk,:) = reshape( myTheta.A(:,:,kk), 1, D*DR );
%         end
%         if K > 1
%             D = pdist( Avec, 'euclidean');        
%             D = squareform(D);
%             D = D(ki,:);
%         else
%             D = 0; % Default distance from self is zero.
%         end
%         D = D./(EPS+max(max(D)) ); 
%         % Ensure D varies from 0 to 1. 
%         %  but is never quite 0 (must be at least EPS)
% end
% 
% D( ~F(jj,:)  ) = Inf;
% 
% ps = 1 - D;
% ps( ~F(jj,:) ) = 0;
% 
% ps = ps./sum( ps );
% 
% if F( jj, ki ) && length( ps( ps>0 ) ) > 1
%     % If object jj possesses feature ki,
%     %    then let prob kj =ki = 2/3
%     %                  kj~=ki = 1/3 (unif. distrib. over all jj's feats)
%     ps( ki ) = 2*( sum( ps( [1:ki-1 ki+1:end] ) ) );
% end
% 
% ps = ps./sum(ps);
% 
% % NOTE: Alternate distance metric
% %sig = eps+min(D(availFeatIDs));
% %ps = exp(-D./sig );

end 