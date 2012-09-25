function algP = defaultMCMCParams_BPHMM()
% Creates a struct encoding the default settings for MCMC inference

algP.Niter = 50;

algP.doSampleFShared = 1;
algP.doSampleFUnique = 1;
algP.doSplitMerge = 0;
algP.doSplitMergeRGS = 0;

algP.SM.featSelectDistr = 'splitBias+margLik';
algP.SM.doSeqUpdateThetaHatOnMerge = 0;

algP.nSMTrials = 5;
algP.nSweepsRGS = 5;

algP.doSampleZ = 1;
algP.doSampleTheta = 1;
algP.doSampleEta   = 1;

algP.BP.doSampleMass = 1;
algP.BP.doSampleConc = 1;
algP.BP.Niter = 10;
algP.BP.var_c = 2;

algP.HMM.doSampleHypers = 1;
algP.HMM.Niter = 20;
algP.HMM.var_alpha = 2;
algP.HMM.var_kappa = 10;

% Reversible Jump Proposal settings
algP.theta.birthPropDistr = 'DataDriven';
algP.theta.minW = 15; 
algP.theta.maxW = 100;

