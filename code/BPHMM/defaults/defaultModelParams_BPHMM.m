function model = defaultModelParams_BPHMM( data )

switch data.getObsType()
    case 'Gaussian'
        model.obsM.precMu   = 1;
        model.obsM.degFree  = 3;
        model.obsM.doEmpCovScalePrior = 0;
        model.obsM.Scoef    = 1;
    case 'AR'
        model.obsM.doEmpCov = 0;
        model.obsM.doEmpCovFirstDiff = 1;
        model.obsM.degFree  = 1;
        if ~isempty( strmatch( '13_29',data.seqNames )  )
            model.obsM.Scoef = 0.5;
        else
            model.obsM.Scoef    = 0.5;
        end
end

% ------------------------------- HMM params
model.hmmM.alpha = 1;
model.hmmM.kappa = 25;
model.hmmM.prior.a_alpha = 0.01;
model.hmmM.prior.b_alpha = 0.01;
model.hmmM.prior.a_kappa = 0.01;
model.hmmM.prior.b_kappa = 0.01;

% ================================================== BETA PROCESS MODEL
% ------------------------------- GAMMA: Mass param for IBP
model.bpM.gamma = 5;
model.bpM.prior.a_mass = 0.01;
model.bpM.prior.b_mass = 0.01;

% ------------------------------- c0   : Concentration param for IBP
model.bpM.c = 1; 
model.bpM.prior.a_conc = 0.01;
model.bpM.prior.b_conc = 0.01;


