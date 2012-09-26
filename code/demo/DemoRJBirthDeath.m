clear all;

doDataDriven = 1;

Ktrue = 4;
if ~exist('data','var')
    [data True] = genToySeqData_Gaussian( Ktrue, 20, 1, 1000 );
end

mP = defaultModelParams_BPHMM(data);
mP.bpM.gamma = 2;

initP = defaultInitMCMC_BPHMM();
algP  = defaultMCMCParams_BPHMM();
outP  = defaultOutputParams_BPHMM({1,1}, 10);

if doDataDriven
    algP.theta.birthPropDistr = 'DataDriven';
else
    algP.theta.birthPropDistr = 'Prior';
end

initP.F.nTotal = 1;

Psi = initBPHMMFresh( data, mP, initP, algP, outP);


for trial = 1:100
    
   [Psi, Stats] = sampleUniqueFeats(Psi, data, algP );
   
   if Stats.ADD.nAccept > 0
      fprintf( 'BIRTH! \n' ); 
      plotGauss( Psi.ThetaM.theta( end ).mu', Psi.ThetaM.theta(end).invSigma, 1 );
   end
   
   Psi = sampleStateSeq( Psi, data );
   Psi.ThetaM = Psi.ThetaM.sampleAllTheta( data, Psi.stateSeq );
   Psi.TransM = Psi.TransM.sampleAllEta( Psi.F, Psi.stateSeq );
   
   if mod(trial,5)==0
       figure(101); clf;
       hold all;
       plotEmissionParams(Psi, data);
       drawnow;
   end
   
end