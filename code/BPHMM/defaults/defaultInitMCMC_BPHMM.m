function initParams = defaultInitMCMC_BPHMM()

% ---------------------------------------------  Initializations
initParams.InitFunc = @initBPHMMFresh;

initParams.F.nUniquePerObj  = 2; 

initParams.z.doPartition    = 1; % 0 := init from prior

