function [Psi, RhoTerms] = ...
    sampleSingleFeatEntry_UniqueRJ( ii, Psi, data, algParams )
% Sample a unique entry in the feature vector of sequence "ii"
% Uses reversible jump to either
%   Create new feature ("birth")
%   Delete cur feature ("death")  
%INPUT
%  ii : integer id of specific sequence to examine
%  Psi : model configuration (includes feat asgns and HMM params)
%  data : SeqObs data object
%           data for sequence ii accessed by call "data.seq(ii)"
%  algParams : param struct specifies details of proposal move
%      algParams.theta.birthPropDistr can be any of:
%      --- 'prior'  : emit param thetaStar drawn from prior
%      --- 'data-driven' : emit param thetaStar draw from data posterior
%OUTPUT
%  Psi : new model config (with potentially new unique features for ii )
%  RhoTerms : some stats about the MH proposal and what kind of move occurs


% ========================================================  UNPACK
F = Psi.F;
TransM = Psi.TransM;
ThetaM = Psi.ThetaM;
gamma  = Psi.bpM.gamma;
c      = Psi.bpM.c;

featureCounts = sum( F, 1 );
availFeatIDs  = find( F(ii,:) > 0  );
K   = size(F,2);
N   = size(F,1);
Kii = length( availFeatIDs );
uniqueFeatIDs = availFeatIDs( featureCounts( availFeatIDs ) == 1  );
uCur = length(uniqueFeatIDs);

% -------------------------------- Eta prop
propEta_ii = TransM.sampleEtaProposal_UniqueBirth( ii );

% -------------------------------  Theta prop
switch algParams.RJ.birthPropDistr
    case {'Prior','prior'}
        choice = 1;
        wstart = 0; wend = 0;
        [thetaStar] = ThetaM.sampleThetaProposal_BirthPrior( );
        propThetaM = ThetaM.insertTheta( thetaStar );
    case {'DataDriven', 'DD', 'datadriven'}
        if isfield( algParams, 'Debug' ) && isfield( algParams.Debug, 'wstart' )        
            wstart = algParams.Debug.wstart;
            wend   = algParams.Debug.wend;
        else
            [wstart, wend, L] = drawRandomSubwindow( data.Ts(ii), algParams.RJ.minW, algParams.RJ.maxW );
        end
        [thetaStar,PPmix, choice] = ThetaM.sampleThetaProposal_BirthDataDriven( ii, data, wstart:wend );
        propThetaM = ThetaM.insertTheta( thetaStar );
end

% Define prob of proposing birth move
%   as deterministic function of the # of unique features
PrBirth = @(uCur)1/2;
qs = buildRJMoveDistr( uCur, Kii, PrBirth );
MoveType = multinomial_single_draw( qs );

if isfield( algParams, 'Debug' )
    MoveType = algParams.Debug.MoveType;
    if MoveType == 0
       if uCur == 0 || Kii==1
           RhoTerms = struct('doBirth',0, 'doAccept',0);
           return; 
       else
           MoveType = 1 + randsample( 1:uCur, 1);
           assert( MoveType <= length(qs), 'Bad choice for debug move selector');
       end
    end
end

if MoveType == 1
    % ----------------------------------------------------  Birth:
    descrStr = 'birth';
    
    kk = K+1;    
    propF_ii = F(ii,:) == 1;
    propF_ii(kk) = 1;
    uNew = uCur + 1;
    
    % ----------------------------------------------- build eta
    % Birth move, keep around *all* of the entries in Pz
    propEta = propEta_ii;
    
    switch algParams.RJ.birthPropDistr
        case {'DataDriven', 'DD', 'datadriven'}   
            if algParams.RJ.doHastingsFactor
            logPrThetaStar_prior = propThetaM.calcLogPrTheta( thetaStar );
            logPrThetaStar_prop  = propThetaM.calcLogPrTheta_MixWithPrior( thetaStar, PPmix );
            logPrTheta_Diff = logPrThetaStar_prior - logPrThetaStar_prop;
            else
                logPrTheta_Diff = -20;
            end
        case {'Prior', 'prior'}
            logPrTheta_Diff = 0;
    end
    
    % Probability of birth in current config
    logQ.moveChoice = log( qs(1) );
    % Probability of killing the last feature  in proposed config
    qsRev = buildRJMoveDistr( uNew, Kii+1, PrBirth );
    logQ_Rev.moveChoice = log( qsRev( end )  );
    
    RhoTerms.activeFeatIDs = kk;
else
    % ----------------------------------------------------  Death:
    descrStr = 'death';
    
    kk =  uniqueFeatIDs( MoveType-1 );
    propF_ii = F(ii,:) == 1;
    propF_ii( kk ) = 0;   
    uNew = uCur - 1;
    
    % ----------------------------------------------- build eta
    jj = find( availFeatIDs == kk );
    keepFeatIDs = [1:jj-1 jj+1:Kii];
    propEta = propEta_ii( keepFeatIDs, keepFeatIDs );
    
    % ----------------------------------------------- build theta
    switch algParams.RJ.birthPropDistr
        case {'DataDriven', 'DD', 'datadriven'}
            thetaStar = propThetaM.theta(kk);
            if algParams.RJ.doHastingsFactor
                logPrThetaKK_prior = propThetaM.calcLogPrTheta( propThetaM.theta(kk) );
                logPrThetaKK_prop  = propThetaM.calcLogPrTheta_MixWithPrior( propThetaM.theta(kk), PPmix );
                logPrTheta_Diff = logPrThetaKK_prop - logPrThetaKK_prior;
            else
                logPrTheta_Diff=-20;
            end
        case {'Prior', 'prior'}
            logPrTheta_Diff = 0;
    end
    
    % Probability of death  in current config
    logQ.moveChoice = log( qs(MoveType) );
    % Probability of birth in proposed config
    qsRev = buildRJMoveDistr( uNew, Kii-1, PrBirth );
    logQ_Rev.moveChoice = log( qsRev( 1 )  );
    
    RhoTerms.activeFeatIDs = kk;
end

% Compute Joint Log Probability of Proposed/Current configurations

% -------------------------------- p( F ) terms
% U ~ Poisson( eta ) where eta = gamma*c/(c + N - 1)
eta = gamma *c/(c + N -1 );
logPrNumFeat_Diff = ( uNew - uCur )*log( eta ) +  gammaln( uCur + 1 ) - gammaln( uNew + 1 );

% -------------------------------- p( x | eta, theta, F)  terms
if isfield( Psi, 'cache' ) && isfield( Psi.cache, 'logSoftEv' )
    logSoftEv = Psi.cache.logSoftEv{ii};
    if strcmp( descrStr, 'birth' )
        logSoftEvStar = propThetaM.calcLogSoftEv( ii, data,  [K+1] );
        logSoftEv( K+1,:) = logSoftEvStar(K+1,:);
    end
else
    logSoftEv = propThetaM.calcLogSoftEv( ii, data, [availFeatIDs K+1] ); 
end
if isfield(Psi,'cache') && isfield( Psi.cache, 'logMargPrObs' )
    logMargPrObs_Cur = Psi.cache.logMargPrObs(ii);
else
    curF_ii = false( size(propF_ii) );
    curF_ii( availFeatIDs ) = true;
    logMargPrObs_Cur = calcLogMargPrObsSeqFAST( logSoftEv( curF_ii, :), propEta_ii( 1:Kii, 1:Kii ) );
end
logMargPrObs_Prop = calcLogMargPrObsSeqFAST( logSoftEv( propF_ii, :), propEta );

% Compute accept-reject ratio:  ( see eq. 15 in BP HMM paper )
log_rho_star = logMargPrObs_Prop - logMargPrObs_Cur  ...
    + logPrNumFeat_Diff ...
    + logPrTheta_Diff ...
    + logQ_Rev.moveChoice - logQ.moveChoice;

RhoTerms.logMargPrObs_Prop = logMargPrObs_Prop;
RhoTerms.logMargPrObs_Cur  = logMargPrObs_Cur;
RhoTerms.logPrThetaDiff = logPrTheta_Diff;
RhoTerms.logQMove.fwd = logQ.moveChoice;
RhoTerms.logQMove.rev = logQ_Rev.moveChoice;
RhoTerms.thetaStar = thetaStar;
RhoTerms.choice = choice;
RhoTerms.window = [wstart wend];

rho = exp(log_rho_star);
assert( ~isnan(rho), 'ERROR: Accept ratio *rho* for unique features should never be NaN')

rho = min( rho, 1 );
doAccept = rand < rho;  % Binary indicator for if cur proposal accepted

RhoTerms.doAccept = doAccept;
RhoTerms.doBirth  = strcmp( descrStr, 'birth' );

if doAccept
    switch descrStr
        case 'birth'
            f_ii_kk = 1;
        case 'death'
            f_ii_kk = 0;
    end
    Psi.F( ii, kk ) = f_ii_kk;
    Psi.ThetaM = propThetaM;
    Psi.TransM = Psi.TransM.setEta( ii, propF_ii, propEta );
    
    if isfield(Psi,'cache') && isfield( Psi.cache, 'logSoftEv' )
        Psi.cache.logMargPrObs(ii) = logMargPrObs_Prop;
        
        if strcmp( descrStr, 'birth' )
            Psi.cache.logSoftEv{ii} = logSoftEv;
        end
    elseif isfield( algParams, 'doAvoidCache' ) && algParams.doAvoidCache
        assert( 1==1 ); % placeholder, skip caching!
    else
        Psi.cache.logMargPrObs(ii) = logMargPrObs_Prop;
        Psi.cache.logSoftEv{ii} = logSoftEv;
    end
    
    if strcmp( descrStr,'death')
      Psi = reallocateFeatIDs( Psi );
    end
end

end % MAIN FUNCTION

% Draw a random subwindow for data-driven proposal
%   First choose window length L
%   Then choose starting position (which can support that length L )
function [wstart, wend, L] = drawRandomSubwindow( TT, minW, maxW )
    Ls = minW:maxW;
    Ls = Ls( Ls >= 1 );
    Ls = Ls( Ls <= TT );
    if isempty(Ls)
        Ls = TT;
    end
    distrLs = ones( size(Ls) );
    L = Ls( multinomial_single_draw( distrLs )  );
    wstart = randi( [1 TT-L+1] );
    wend   = wstart + L - 1;
end
