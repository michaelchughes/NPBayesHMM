function [Psi, algParams, outParams] = initBPHMMCheat( data, model, initParams, algParams, outParams )
% Initialize BP-HMM hidden variables for MCMC sampling as follows
%   depends on parameter settings.Cheat.nRepeats
%          which will create "nRepeats" copies of each "true" feature
%             by dividing data items into nRepeats blocks
%                  and giving each block its own set of features
%   3) sample TS from prior given F, z
%   4) sample emit params Theta from posterior given stateSeq, data


% ------------------------------------------------ Build cheat config F, z
if exist('initParams','var') && isfield( initParams, 'Cheat') && isfield( initParams.Cheat, 'nRepeats')
    nRepeats = initParams.Cheat.nRepeats;
else
    nRepeats = 1;
end
blockSizes = getEqualSizePartitions( data.N, nRepeats );
cSizes = cumsum( blockSizes );
baseFeatID = 0;
for bb = 1:length( blockSizes )
   if bb == 1
        objIDs = 1:blockSizes(1);
   else
        objIDs = cSizes(bb-1)+1:cSizes(bb);
   end
   
   uBlock = [];
   for ii = objIDs
       uBlock = [uBlock data.zTrue(ii) ];
   end
   uFeatIDs = unique( uBlock );
      
   for ii = objIDs
       for uu = 1:length( uFeatIDs )
           ttINDS = data.zTrue(ii) == uFeatIDs(uu);
           if sum( ttINDS ) > 0
                F( ii, baseFeatID+uu ) = 1;
                stateSeq(ii).z( ttINDS ) = baseFeatID + uu;
           end
       end
   end
   baseFeatID = baseFeatID + length(uFeatIDs);
end

% ------------------------------------------------- Init Trans Struct (Eta)
% Sample the transition weights eta and init. state weights eta_init
%   from the priors on these distributions
TransM = HMMTransModel(  size(F,1), size(F,2) );
TransM = TransM.setPrior( model.hmmM.alpha, model.hmmM.kappa,...
        model.hmmM.prior.a_alpha, model.hmmM.prior.b_alpha, ...
        model.hmmM.prior.a_kappa, model.hmmM.prior.b_kappa );
TransM = TransM.sampleAllEta( F, stateSeq );

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
