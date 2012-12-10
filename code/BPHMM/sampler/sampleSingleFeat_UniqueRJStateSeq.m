function [Psi, RhoTerms] = ...
    sampleSingleFeat_UniqueRJStateSeq( ii, Psi, data, algParams )
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
%OUTPUT
%  Psi : new model config (with potentially new unique features for ii )
%  RhoTerms : some stats about the MH proposal and what kind of move occurs

doDebug=0;

% ========================================================  UNPACK
F = Psi.F;

propStateSeq = Psi.stateSeq;
propF = F;

gamma  = Psi.bpM.gamma;
c      = Psi.bpM.c;

featureCounts = sum( F, 1 );
availFeatIDs  = find( F(ii,:) > 0  );
K   = size(F,2);
N   = size(F,1);
Kii = length( availFeatIDs );
uniqueFeatIDs = availFeatIDs( featureCounts( availFeatIDs ) == 1  );
uCur = length(uniqueFeatIDs);

% -------------------------------- DETERMINISTIC Eta prop
propTransM = Psi.TransM.getAllEta_PriorMean( ii, [F zeros(N,1)], K+1 );
EtaHatAll = propTransM.seq(ii);

% -------------------------------  DETERMINISTIC Theta prop
ThetaHat = Psi.ThetaM;
for kk = availFeatIDs
    PN = ThetaHat.getPosteriorParams( ThetaHat.Xstats(kk) );
    ThetaHat.theta(kk) = Psi.ThetaM.getTheta_Mean( PN );
end

switch algParams.RJ.birthPropDistr
    case {'Prior','prior'}
        wstart = 0; wend = 0;
        [thetaStar] = ThetaHat.getTheta_Mean( );
        ThetaHat = ThetaHat.insertTheta( thetaStar );
    case {'DataDriven', 'DD', 'datadriven'}
        [wstart, wend, L] = drawRandomSubwindow( data.Ts(ii), algParams.RJ.minW, algParams.RJ.maxW );
        if strcmp( class(data), 'ARSeqData' )
            X = data.seq(ii);
            Xprev = data.prev(ii);
            PN = ThetaHat.getPosteriorParams( ThetaHat.getXSuffStats( X(:,wstart:wend), Xprev(:,wstart:wend) ) );
        else
            X = data.seq(ii);
            PN = ThetaHat.getPosteriorParams( ThetaHat.getXSuffStats( X(:,wstart:wend) ) );
        end
        thetaStar = ThetaHat.getTheta_Mean(PN);
        ThetaHat = ThetaHat.insertTheta( thetaStar );
end

seqSoftEv = ThetaHat.calcLogSoftEv(  ii, data, [availFeatIDs K+1] );

% Define prob of proposing birth move
%   as deterministic function of the # of unique features
PrBirth = @(uCur)1/2;
qs = buildRJMoveDistr( uCur, Kii, PrBirth );
MoveType = multinomial_single_draw( qs );

if isfield( algParams, 'Debug' )
    MoveType = algParams.Debug.MoveType;
end

if MoveType == 1
    % ----------------------------------------------------  Birth:
    descrStr = 'birth';
    
    kk = K+1;    
    propF_ii = F(ii,:) == 1;
    propF_ii(kk) = 1;
    propF(ii,kk)=1;
    uNew = uCur + 1;
    
    % ----------------------------------------------- build eta
    % Birth move, keep around *all* of the entries in Pz
    EtaHat = EtaHatAll;
    
    % ----------------------------------------------- sample proposed z_ii
    propFeatIDs = [availFeatIDs K+1];
    [propStateSeq(ii).z, logQ.z] = sampleSingleStateSeq_WithSoftEv( ii, EtaHat, seqSoftEv(propFeatIDs,:) );
    
    propThetaHat = ThetaHat.decXStats( ii, data, Psi.stateSeq, availFeatIDs );
    propThetaHat = propThetaHat.incXStats( ii, data, propStateSeq, propFeatIDs );
    for jj = propFeatIDs
        PN = propThetaHat.getPosteriorParams( propThetaHat.Xstats(jj) );
        propThetaHat.theta(jj) = propThetaHat.getTheta_Mean( PN );
    end
    seqSoftEvRev = propThetaHat.calcLogSoftEv(  ii, data, [availFeatIDs] );
    seqSoftEvRev = seqSoftEvRev(availFeatIDs,:);
    % ----------------------------------------------- reverse to original z
    EtaHatOrig.availFeatIDs = availFeatIDs;
    EtaHatOrig.eta = EtaHatAll.eta( 1:Kii, 1:Kii);
    [~, logQ_Rev.z] = sampleSingleStateSeq_WithSoftEv( ii, EtaHatOrig, seqSoftEvRev, Psi );
    
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
    propF(ii,kk)=0;
    uNew = uCur - 1;
    
    % ----------------------------------------------- build eta
    jj = find( availFeatIDs == kk );
    keepFeatIDs = [1:jj-1 jj+1:Kii];
    EtaHat.availFeatIDs = availFeatIDs(keepFeatIDs);
    EtaHat.eta = EtaHatAll.eta( keepFeatIDs, keepFeatIDs );
    
    % ----------------------------------------------- sample proposed z_ii
    [propStateSeq(ii).z, logQ.z] = sampleSingleStateSeq_WithSoftEv( ii, EtaHat, seqSoftEv( availFeatIDs(keepFeatIDs),:) );
   
    propThetaHat = ThetaHat.decXStats( ii, data, Psi.stateSeq, availFeatIDs );
    propThetaHat = propThetaHat.incXStats( ii, data, propStateSeq, availFeatIDs );
    for jj = availFeatIDs(keepFeatIDs)
        PN = propThetaHat.getPosteriorParams( propThetaHat.Xstats(jj) );
        propThetaHat.theta(jj) = propThetaHat.getTheta_Mean( PN );
    end
    seqSoftEvRev = propThetaHat.calcLogSoftEv(  ii, data, [availFeatIDs(keepFeatIDs) K+1] );
    seqSoftEvRev = seqSoftEvRev( [availFeatIDs(keepFeatIDs) K+1],:);
    
    % ----------------------------------------------- reverse to original z
    EtaHatOrig.availFeatIDs = availFeatIDs;
    EtaHatOrig.eta = EtaHatAll.eta( 1:Kii, 1:Kii);
    [~, logQ_Rev.z] = sampleSingleStateSeq_WithSoftEv( ii, EtaHatOrig, seqSoftEvRev, Psi );
    
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

% -------------------------------- p( z_ii | F ) term
logPrZ_Prop = Psi.TransM.calcMargPrStateSeq( propF, propStateSeq, ii );
logPrZ_Cur  = Psi.TransM.calcMargPrStateSeq( F, Psi.stateSeq, ii );

% -------------------------------- p( x | z, F)  terms
if MoveType==1
    logPrObs_Prop = propThetaHat.calcMargPrData( data, propStateSeq, propFeatIDs );
else
    logPrObs_Prop = propThetaHat.calcMargPrData( data, propStateSeq, availFeatIDs(keepFeatIDs) );
end
% Only concerned with the active features!
logPrObs_Cur  = Psi.ThetaM.calcMargPrData([], [], availFeatIDs);

logQHastings = logQ_Rev.z - logQ.z ...
    + logQ_Rev.moveChoice - logQ.moveChoice;

if algParams.doAnneal
    if Psi.invTemp == 0 && isinf(logQHastings)
        logQHastings = -Inf; % always want to ignore in this case!
        % this is a sign of something seriously bad with construction
    else
        logQHastings = Psi.invTemp * logQHastings;
    end
end

% Compute accept-reject ratio:  ( see eq. 15 in BP HMM paper )
log_rho_star = logPrObs_Prop - logPrObs_Cur ...
    + logPrNumFeat_Diff ...
    + logPrZ_Prop - logPrZ_Cur ...
    + logQHastings;

RhoTerms.thetaStar = thetaStar;
RhoTerms.window = [wstart wend];

rho = exp(log_rho_star);
assert( ~isnan(rho), 'ERROR: Accept ratio *rho* for unique features should never be NaN')

rho = min( rho, 1 );
doAccept = rand < rho;  % Binary indicator for if cur proposal accepted

RhoTerms.doAccept = doAccept;
RhoTerms.doBirth  = strcmp( descrStr, 'birth' );

% if doDebug
%     propPsi.F    = propF;
%     propPsi.stateSeq = propStateSeq;
%     propPsi.TransM   = Psi.TransM;
%     propPsi.TransM.seq(ii) = EtaHat;
%     propPsi.ThetaM   = propThetaHat;
%     propPsi.bpM  = Psi.bpM;
%     
%     lPP = calcJointLogPr_BPHMMState( propPsi, data);
%     lPC = calcJointLogPr_BPHMMState( Psi,     data);
%     assert( allEq( lPP.obs-lPC.obs, logPrObs_Prop-logPrObs_Cur), 'bad' );
%     assert( allEq( lPP.z -lPC.z , logPrZ_Prop-logPrZ_Cur), 'bad' );
% end

if doAccept
    switch descrStr
        case 'birth'
            f_ii_kk = 1;
        case 'death'
            f_ii_kk = 0;
    end
    Psi.F( ii, kk ) = f_ii_kk;
    Psi.stateSeq = propStateSeq;
    Psi.ThetaM = propThetaHat;
    Psi.TransM = Psi.TransM.setEta( ii, propF_ii, EtaHat.eta );
    
    if isfield( Psi, 'cache')
        Psi = rmfield( Psi, 'cache');
    end
    %if isfield(Psi,'cache') && isfield( Psi.cache, 'logSoftEv' )        
    %    if strcmp( descrStr, 'birth' )
    %        Psi.cache.logSoftEv{ii} = seqSoftEv;
    %    end
    %else
    %    Psi.cache.logSoftEv{ii} = seqSoftEv;
    %end
    
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
