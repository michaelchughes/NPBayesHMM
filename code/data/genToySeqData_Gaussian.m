function [data, PsiTrue] = genToySeqData_Gaussian( nStates, nDim, N, T, pIncludeFeature)
% INPUTS ----------------------------------------------------------
%    nStates = # of available Markov states
%    nDim = number of observations at each time instant
%    N = number of time series objects
%    T = length of each time series
% OUTPUT ----------------------------------------------------------
%    data  :  SeqData object

% ------------------------------- Remember old state to use again afterward
curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to default value with SEED 0
%       so that we always get same synth data regardless of when called
reset( RandStream.getGlobalStream(), 0);

if T < 0
    doVaryLength = 1;
    T = abs(T);
else
    doVaryLength = 0;
end

if ~exist( 'pIncludeFeature', 'var' )
    pIncludeFeature = 0.75;
end
pSelfTrans = 1-(2*nStates)/T;

% Create initial state distribution (uniform)
Pz_init = ones(1, nStates);

% Create state-to-state transition matrix Pz
Pz = zeros(  nStates, nStates );
for k = 1:nStates
    Pz(k,k) = pSelfTrans;
    Pz(k, [1:k-1 k+1:end] ) = (1-pSelfTrans)/(nStates-1);
end
doRowsSumToOne = ( sum(Pz,2) - ones(size(Pz,1),1) ) <= 1e-10;
assert( all(doRowsSumToOne), 'ERROR: Not a valid transition distr.' );

% Create state-specific emission params Px
%   Means are evenly spaced around the unit circle
%   Covariances are aligned so major diagonal of ellipsoid
%        points toward the origin
Px.Mu = zeros( nStates, nDim );
Px.Sigma = zeros( nDim, nDim, nStates );
V = 0.05;
ts = linspace( -pi, pi, nStates+1 );
xs = cos(ts);
ys = sin(ts);

S(:,:,1) = [V 0; 0 .1*V];
S(:,:,2) = [0.5*V .4*V; .4*V .5*V];
S(:,:,3) = [0.1*V 0; 0 V];
S(:,:,4) = [0.5*V -.4*V; -.4*V .5*V];
if nDim == 1
    Px.Mu = linspace( -1*V, 1*V, nStates )';
    for kk = 1:nStates
        Px.Sigma(:,:,kk) = V;
    end
elseif nDim == 2
    Px.Mu = [xs(1:end-1)' ys(1:end-1)'];
    
    if nStates == 8
        Px.Sigma(:,:,1:4) = S(:,:,1:4);
        Px.Sigma(:,:,5:8) = S(:,:,1:4);
    elseif nStates == 4
        Px.Sigma(:,:,1:4) = S(:,:,[1 3 1 3] );
    else
        for kk = 1:nStates
            Px.Sigma(:,:, mod(kk,4)+1 ) = S(:,:,mod(kk,4)+1);
        end
    end
else
    Px.Mu = [xs(1:end-1)' ys(1:end-1)'  zeros( nStates, nDim-2 ) ];
    
    if nDim > 10
       nExtras = floor( (nDim-1) / 10 ); 
       Px.Mu(:, 10*(1:nExtras) ) = repmat( 10*sqrt(V)*xs(1:end-1)', 1, nExtras );
       Px.Mu(:, 10*(1:nExtras)+1 ) = repmat( 10*sqrt(V)*ys(1:end-1)', 1, nExtras );
    end
    
    for kk = 1:nStates
       Px.Sigma(1:2, 1:2, kk) = S(:, :, mod(kk-1, 4)+1  ); 
       Px.Sigma(3:end, 3:end,kk) = V*eye( nDim-2 );
    end
end


% Build time series
data = SeqData();
F = zeros( N, nStates );
for i = 1:N   
    if doVaryLength
        Ti = poissrnd(T);
    else
        Ti = T;
    end
      
    % Draw subset of states that this time-series exhibits
    mask = rand( 1, nStates ) < pIncludeFeature;
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( nStates, 1);
        mask(  kk  ) = 1;
    end
    if i == 1
        mask=true(1,nStates);
    end
    
    F(i,:) = mask;
    
    zTrue = zeros(1,Ti);
    X = zeros( nDim, Ti );
    for t = 1:Ti
        % ---------------------------------------------- Assign true label
        if t == 1
            zcur = multinomial_single_draw( mask.*Pz_init );
        else
            zprev = zTrue(t-1);
            zcur = multinomial_single_draw( mask.*Pz(zprev,:) );
        end
        zTrue(t) = zcur;

        % ---------------------------------------- Assign emissions
       X(:,t) = mvnrnd( Px.Mu(zcur,:),  Px.Sigma(:,:,zcur) );
           
    end
    data = data.addSeq( X, num2str(i), zTrue );
end

% ---------------------------------------------------------  Reset stream
curStream = RandStream.getGlobalStream();
curStream.State = entryState;

PsiTrue.F = zeros(N, nStates);
for ii = 1:N
    PsiTrue.F(ii, unique( data.zTrue(ii) ) ) = 1;
end
for kk = 1:nStates
    PsiTrue.theta(kk).mu = Px.Mu(kk,:);
    PsiTrue.theta(kk).invSigma = inv( Px.Sigma(:,:,kk) );
end
PsiTrue.Pz = Pz;
PsiTrue.z = zTrue;

end % main function



