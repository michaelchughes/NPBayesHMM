function ps = getSplitMergeFeatChoiceProbs( F, stateSeq, data_struct, model, jj, ki )
%function ps = getSplitMergeFeatChoiceProbs( F, jj, ki )
ps = F(jj,:);

availFeatIDs = find( F(jj,:) );
relFeatIDs = union( availFeatIDs, ki);

%Ustats  = getStateSeqSuffStats_C1( stateSeq, F, data_struct, model, 1:size(F,1), 1:size(F,2) );
%Ustats.Nkv = Ustats.Nkv(relFeatIDs,:);
%myTheta = getTheta_PosteriorMean( myTheta, Ustats, model.obsModel, relFeatIDs );

Xstats = getXSuffStats( F, stateSeq, data_struct, model, 1:size(F,1),  relFeatIDs );
myTheta = setThetaToPosteriorMean(  [], Xstats, model.obsModel, relFeatIDs );
        

EPS = 1e-9;
switch model.obsModel.type
    case 'Multinomial'
        pHat = exp( myTheta.logp )';
        pHat = bsxfun( @rdivide, pHat, sum(pHat,1) );
        D = calcChiSqDistanceBetter( pHat(:,ki), pHat );
    case 'Gaussian'
        if size( myTheta.Mu, 1) > 1
            D = pdist( myTheta.Mu, 'euclidean');        
            D = squareform(D);
            D = D(ki,:);
        else
            D = 0; % Default distance from self is zero.
        end
        D = D./(EPS+max(max(D)) ); % Ensure D varies from 0 to 1. 
        
    case 'AR-Gaussian'
        [D DR K] = size( myTheta.A );
        
        %Avec = -100*ones(K,D);
        Avec = zeros(K, D*DR);
        for kk = relFeatIDs
            %Sig = myTheta.invSigma(:,:,kk) \ eye(D);
            %Avec(kk,:) = diag( Sig );
            Avec(kk,:) = reshape( myTheta.A(:,:,kk), 1, D*DR );
        end
        if K > 1
            D = pdist( Avec, 'euclidean');        
            D = squareform(D);
            D = D(ki,:);
        else
            D = 0; % Default distance from self is zero.
        end
        D = D./(EPS+max(max(D)) ); 
        % Ensure D varies from 0 to 1. 
        %  but is never quite 0 (must be at least EPS)
end

D( ~F(jj,:)  ) = Inf;

ps = 1 - D;
ps( ~F(jj,:) ) = 0;

ps = ps./sum( ps );

if F( jj, ki ) && length( ps( ps>0 ) ) > 1
    % If object jj possesses feature ki,
    %    then let prob kj =ki = 2/3
    %                  kj~=ki = 1/3 (unif. distrib. over all jj's feats)
    ps( ki ) = 2*( sum( ps( [1:ki-1 ki+1:end] ) ) );
end

ps = ps./sum(ps);

% NOTE: Alternate distance metric
%sig = eps+min(D(availFeatIDs));
%ps = exp(-D./sig );

end 