function Psi = unpackBPHMMState( sPsi, data, model );
% unpackBPHMMState : convert stored sampler state into fully usable one
% EXPLANATION:
% HMM parameters for transitions + emissions are represented as objects
% but this takes up lots of room when storing results for many iterations
%   so we store the objects as raw numeric parameters
%   and then rebuild by calling this function

Psi.F = sPsi.F;
Psi.stateSeq = sPsi.stateSeq;
Psi.bpM = model.bpM;
Psi.bpM.gamma = sPsi.gamma;
Psi.bpM.c     = sPsi.c;

F = Psi.F;
TransM = HMMTransModel(  size(F,1), size(F,2) );

if isfield( sPsi, 'alpha' )
    alph = sPsi.alpha;
    kapp = sPsi.kappa;
else
    alph = model.hmmM.alpha;
    kapp = model.hmmM.kappa;
end

TransM = TransM.setPrior( alph, kapp, ...
    model.hmmM.prior.a_alpha, model.hmmM.prior.b_alpha, ...
        model.hmmM.prior.a_kappa, model.hmmM.prior.b_kappa );
for ii = 1:size( F, 1)
    TransM = TransM.setEta( ii, F(ii,:), sPsi.Eta(ii).eta );
end

switch data.getObsType()
    case 'Gaussian'
      ThetaM = SeqObsModel_Gaussian( 0, data.D );
    case 'AR'
      ThetaM = SeqObsModel_ARGaussian( 0, data.D, data.R );
end
ThetaM = ThetaM.setPrior( data , model.obsM );
for kk = 1: size( Psi.F,2 )
    ThetaM = ThetaM.insertTheta( sPsi.theta(kk) );
end

Psi.TransM = TransM;
Psi.ThetaM = ThetaM;

end