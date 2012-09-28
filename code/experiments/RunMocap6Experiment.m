function [] = RunMocap6Experiment( jobID, taskID, infName, initName, TimeLimit )
%INPUT
%  jobID : integer/name of job
%  taskID : integer/name of task
%  infName  : {'Prior','DD','SM'} indicates type of sampler     
%  initName : {'one', 'unique5'}  indicates initialization
%  TimeLimit : # of seconds to allow for MCMC


dataP = {'Mocap6'};
modelP = {};

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
    case 'SMnoqrev'               
        algP = {'doSampleFUnique', 0, 'doSMNoQRev', 1}; 
    case 'SMnoqrev+DD'               
        algP = {'doSampleFUnique', 1, 'doSMNoQRev', 1, 'theta.birthPropDistr', 'DataDriven'}; 
    otherwise
        error( 'unrecognized inference type' );
end
algP(end+1:end+2) = {'TimeLimit', TimeLimit};

switch initName
    case 'cheat'
        initP = {'InitFunc', @initBPHMMCheat, 'Cheat.nRepeats', 1};
    case 'one'
        initP = {'F.nTotal', 1};
    case 'unique5'
        initP = {'F.nUniquePerObj', 5};
end

runBPHMM( dataP, modelP, {jobID, taskID}, algP, initP );
