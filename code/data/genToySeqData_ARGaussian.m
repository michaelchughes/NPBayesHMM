function [data] = genToySeqData_ARGaussian( nStates, nDim, N, T, R)
% INPUTS ----------------------------------------------------------
%    nStates = # of available Markov states
%    nDim = number of observations at each time instant
%    N = number of time series objects
%    T = length of each time series
%    R = order of the autoregressive process
% OUTPUT ----------------------------------------------------------
%    data  :  ARSeqData object
%               note that each sequence will *actually* have T-R
%               observations, since we need "R" to properly
%                define the likelihood of the first "kept" observation

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

pIncludeFeature = 0.75;
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

as = linspace( -0.9, 0.9, nStates );
A  = zeros( nDim, nDim*R, nStates);

nu = nDim+2;   % inverse Wishart degrees of freedom
meanSigma = 0.5*eye(nDim); % inverse Wishart mean covariance matrix
nu_delta = (nu-nDim-1)*meanSigma;


My2DSigs(:,:,1) = [15 0; 0 0.5];
My2DSigs(:,:,2) = [10 5; 5 3];
My2DSigs(:,:,3) = [0.5 0; 0 15];
My2DSigs(:,:,4) = [3 -5; -5 10];

MyVarCoefs = repmat( [0.001 0.1 1 2 3 4 5 10 15], 1, ceil(nDim/5) );

if R > 1
muteVec = logspace( 0, -R/2, R);
else
muteVec = 1;
end
for kk = 1:nStates
    
    for dd = 1:nDim
         A(dd, [dd:nDim:nDim*R], kk ) = as(kk) .* muteVec;
    end
    
    if nDim == 1
        % IW(0.5, 3)
        Sigma(:,:,kk) = iwishrnd( nu_delta,  nu  );
    else
        Sigma(1:2, 1:2, kk) = My2DSigs(:,:, mod(kk-1,4)+1 );
        Sigma(3:nDim, 3:nDim, kk) = diag( randsample( MyVarCoefs, nDim-2 ) );
    end
end

Px.A = A;
Px.Sigma = Sigma;

% Build time series
data = ARSeqData( R );
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
    F(i,:) = mask;
    
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( nStates, 1);
        %kk = multinomial_single_draw( ones( Ktrue,1 ) );
        mask(  kk  ) = 1;
    end
    F(i,:) = mask;
    
    
    zTrue = zeros(1,Ti);
    X = zeros( nDim, Ti );
    xprev = zeros( nDim*R, 1 );
    for t = 1:Ti
        % ---------------------------------------------- Assign true label
        if t == 1
            zcur = multinomial_single_draw( mask.*Pz_init );
        else
            zprev = zTrue(t-1);
            zcur = multinomial_single_draw( mask.*Pz(zprev,:) );
        end
        zTrue(t) = zcur;
        
        % ---------------------------------------- Assign emitted data
 
        X(:,t) = mvnrnd( Px.A(:,:,zcur)*xprev,  Px.Sigma(:,:,zcur) );
        xprev = [X(:,t); xprev( 1:(end-nDim) ) ];
        
    end
    data = data.addSeq( X, num2str(i), zTrue );
end

% ---------------------------------------------------------  Reset stream
curStream = RandStream.getGlobalStream();
curStream.State = entryState;

end % main function



