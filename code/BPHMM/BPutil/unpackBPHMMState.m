function Psi = unpackBPHMMState( sPsi, data, model );

Psi.F = sPsi.F;
Psi.stateSeq = sPsi.stateSeq;
Psi.bpM = model.bpM;
Psi.bpM.gamma = sPsi.gamma;
Psi.bpM.c     = sPsi.c;

F = Psi.F;
TransM = HMMTransModel(  size(F,1), size(F,2) );

TransM = TransM.setPrior( model.hmmM.alpha, model.hmmM.kappa,...
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