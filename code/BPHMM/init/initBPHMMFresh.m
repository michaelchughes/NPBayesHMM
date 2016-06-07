function [Psi, algParams, outParams] = initBPHMMFresh( data, model, initParams, algParams, outParams )
% Initialize BP-HMM hidden variables for MCMC sampling as follows
%   1) fix F to user-specified configuration
%          each object has U unique features
%   2) sample TS from prior given F
%   3) sample state asgns Z given F, TS
%             object i has U_i features in F
%             divide z_i into U_i equal sized contiguous blocks
%                 assign each to a distinct element of the U_i features
%   4) sample emit params Theta from posterior given stateSeq, data


% ------------------------------- Remember old state to use again afterward
curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to associate with specific *task*
reset( RandStream.getGlobalStream(), outParams.taskID );

if outParams.doPrintHeaderInfo
    fprintf( 'Psi MCMC Chain State:\n' );
end

% --------------------------------------------------- Init Feature Matrix
nObj = data.N;
if isfield( initParams.F, 'nTotal' )
    K =  initParams.F.nTotal;
    F = ones( nObj, K );    
    if outParams.doPrintHeaderInfo        
        fprintf( '\t F : %d global features shared by all objects \n', K );
    end
elseif isfield( initParams.F, 'nUniquePerObj' )
    Kii = initParams.F.nUniquePerObj;
    F = zeros( nObj, Kii * nObj );
    for ii = 1:nObj
        F(ii, (ii-1)*Kii + [1:Kii] ) = 1;
    end
    if outParams.doPrintHeaderInfo        
        fprintf( '\t F : each obj. assigned %d unique features \n', Kii );
    end    
end
% ------------------------------------------------- Init Trans Struct (Eta)
% Sample the transition weights eta and init. state weights eta_init
%   from the priors on these distributions
TransM = HMMTransModel(  size(F,1), size(F,2) );
TransM = TransM.setPrior( model.hmmM.alpha, model.hmmM.kappa,...
    model.hmmM.prior.a_alpha, model.hmmM.prior.b_alpha, ...
        model.hmmM.prior.a_kappa, model.hmmM.prior.b_kappa );
TransM = TransM.sampleAllEta( F );
if outParams.doPrintHeaderInfo
    fprintf( '\t Eta : sampled from prior given avail. states for each obj. \n' );
end

% ---------------------------------------------------   Init stateSeq
stateSeq = repmat( struct('z', zeros(1,100) ),  1, size(F,1) );
if initParams.z.doPartition
    if outParams.doPrintHeaderInfo
        fprintf( '\t StateSeq : fixed to sequence of contiguous blocks \n' );
    end
    for ii = 1:data.N
        T = data.Ts(ii);
        stateSeq(ii).z = zeros( 1, T, 'uint16' );
        availFeatIDs = TransM.seq(ii).availFeatIDs;
        K_ii = length( availFeatIDs );
        blockSizes = getEqualSizePartitions( T, K_ii );
        bstart = 0;
        for kk = 1:K_ii
            blockIDs = bstart+1:bstart+blockSizes(kk);
            stateSeq(ii).z( blockIDs ) = availFeatIDs( kk );
            bstart = bstart + blockSizes(kk);
        end
    end
else
    fprintf( '\t StateSeq : sampled from prior given Eta \n' );
    for ii = 1:data.N
        T = data.Ts(ii);
        stateSeq(ii).z = zeros( 1, T, 'uint16' );
        availFeatIDs = TransM.seq(ii).availFeatIDs;
        pi_init = TransM.seq(ii).pi_init;
        pi_z    = TransM.seq(ii).pi_z;
        for t=1:T
            if (t == 1)
                Pz = pi_init;
            else
                jprev = jcur;
                Pz = pi_z( jprev,: );
            end
            jcur = multinomial_single_draw( Pz );
            stateSeq(ii).z(t) = availFeatIDs( jcur );            
        end
    end
end

% -------------------------------------------------  Init Emission theta
switch data.getObsType()
    case 'Gaussian'
      ThetaM = SeqObsModel_Gaussian( size(F,2), data.D ); 
    case 'AR'
      ThetaM = SeqObsModel_ARGaussian( size(F,2), data.D, data.R);
end

ThetaM = ThetaM.setPrior( data , model.obsM );
ThetaM = ThetaM.sampleAllTheta( data, stateSeq );

% ---------------------------- Repack
Psi.F = F;
Psi.TransM = TransM;
Psi.ThetaM = ThetaM;
Psi.stateSeq = stateSeq;
Psi.bpM = model.bpM;


% ---------------------------------------------------------  Reset stream
curStream = RandStream.getGlobalStream();
curStream.State = entryState;
