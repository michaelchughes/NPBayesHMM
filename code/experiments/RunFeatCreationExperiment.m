function [] = RunFeatCreationExperiment( jobID, taskID, dataName, infName, TimeLimit, doFixHypers )
%INPUT
%  jobID : integer/name of job
%  taskID : integer/name of task
%  dataName : {Gaussian,AR} indicates emission type of toy data
%  infName  : {'Prior','DD','SM'} indicates type of sampler     
%  TimeLimit : # of seconds to allow for MCMC

if ~exist( 'doFixHypers', 'var' )
    doFixHypers = 'BP+HMM'; 
end

switch dataName
    case {'Multinomial', 'Mult'}
        dataP = {'Synth', 'nObj', 100, 'nStates', 8, 'V', 1000, 'obsDim', -5, 'T', 200, 'pEmitFavor', 0.9};
    case {'Gaussian'}        
        dataP = {'SynthGaussian', 'nObj', 100, 'nStates', 8, 'obsDim', 8, 'T', 200 };
        modelP = {'bpM.gamma', 2};
    case {'AR', 'AR-Gaussian'}        
        dataP = {'SynthAR', 'nObj', 100, 'nStates', 8, 'obsDim', 5, 'T', 200, 'R', 1};
        modelP = {'bpM.gamma', 2};
end

initP = {'F.nTotal', 1};

TimeLimit = force2double( TimeLimit );
switch infName
    case 'Prior'
        algP = {'doSampleFUnique', 1, 'doSplitMerge', 0, 'theta.birthPropDistr', 'Prior'};         
    case {'DD', 'DataDriven'}
        algP = {'doSampleFUnique', 1, 'doSplitMerge', 0, 'theta.birthPropDistr', 'DataDriven'};            
    case 'SM'                
        algP = {'doSampleFUnique', 0, 'doSplitMerge', 1};         
    case 'SM+DD'               
        algP = {'doSampleFUnique', 1, 'doSplitMerge', 1, 'theta.birthPropDistr', 'DataDriven'};            
    otherwise
        error( 'unrecognized inference type' );
end

algP(end+1:end+2) = {'TimeLimit', TimeLimit};

if strcmp( doFixHypers, 'c' )
    algP(end+1:end+2) = {'BP.doSampleConc', 0};
elseif strcmp( doFixHypers, 'BP' ) 
    algP(end+1:end+4) = {'BP.doSampleConc', 0, 'BP.doSampleMass', 0};
elseif strcmp( doFixHypers, 'BP+HMM' ) 
    algP(end+1:end+6) = {'BP.doSampleConc', 0, 'BP.doSampleMass', 0, ...
                         'HMM.doSampleHypers', 0};
end

if exist( 'saveEvery', 'var' )
    outP = {jobID, taskID, 'saveEvery', saveEvery, 'printEvery', 10};
else
    outP = {jobID, taskID};
end

runBPHMM( dataP, modelP, outP, algP, initP );
