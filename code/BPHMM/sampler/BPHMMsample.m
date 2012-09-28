function [Psi, Stats] = BPHMMsample( Psi, data, algP)
% Perform ONE iteration of sampling on the BPHMM model,
%  aggregating stats about sampler performance (accept rates) along the way
% Sweeps through the following MCMC moves:
%    shared features F     (MH updates, marg out z)
%    unique features F     (RJ birth/death, marg out z)
%    state sequences z     (Gibbs block sampling via dyn programming)
%    split merge     F,z   (sequential allocation proposals, MH accepted)
%    HMM trans.      eta   (Gibbs updates given z)
%    HMM emission    theta (Gibbs updates given z
%    HMM hypers      alph/kappa  (MH proposals via Gamma random walk)
%    BP  hypers      gamma/c     (MH proposals / Gibbs updates)

if algP.doSampleFShared
    [Psi, Stats.FMH] = sampleSharedFeats( Psi, data );
end

if algP.doSampleFUnique
    [Psi, Stats.FRJ] = sampleUniqueFeats( Psi, data, algP );
end

if algP.doSampleZ
    Psi = sampleStateSeq( Psi, data );
    Psi.ThetaM = Psi.ThetaM.updateAllXSuffStats( horzcat(Psi.stateSeq(:).z), data );
end

if algP.doSplitMerge
    SM.ADD.nAccept=0;
    SM.ADD.nTotal =0;
    SM.DEL.nAccept=0;
    SM.DEL.nTotal =0;
    for trial = 1:algP.nSMTrials
        [nPsi, tS] = sampleSplitMerge_SeqAlloc( Psi, data, algP );
        Psi=nPsi;
        SM.(tS.moveDescr).nAccept = SM.(tS.moveDescr).nAccept + tS.nAccept;
        SM.(tS.moveDescr).nTotal  = SM.(tS.moveDescr).nTotal + 1;
    end
    Stats.SM = SM;
elseif algP.doSMNoQRev
    SM.ADD.nAccept=0;
    SM.ADD.nTotal =0;
    SM.DEL.nAccept=0;
    SM.DEL.nTotal =0;
    for trial = 1:algP.nSMTrials
        [nPsi, tS] = sampleSplitMerge_NoQRev( Psi, data, algP );
        Psi=nPsi;
        SM.(tS.moveDescr).nAccept = SM.(tS.moveDescr).nAccept + tS.nAccept;
        SM.(tS.moveDescr).nTotal  = SM.(tS.moveDescr).nTotal + 1;
    end
    Stats.SM = SM;
end

if algP.doSampleEta
    Psi.TransM = Psi.TransM.sampleAllEta( Psi.F, Psi.stateSeq );
end

if algP.doSampleTheta
    Psi.ThetaM = Psi.ThetaM.sampleAllTheta( data, Psi.stateSeq );
end

if algP.HMM.doSampleHypers
    [Psi, Stats.HMMalpha, Stats.HMMkappa] = sampleHMMhypers( Psi, algP );
end

if algP.BP.doSampleMass || algP.BP.doSampleConc
    [Psi, Stats.BPconc] = sampleIBPhypers(Psi, algP);
end

end % main function