function [ChainHist] = resumeBPHMM( outParams, algParams, initParams )
% runBPHMM
% User-facing entry function for configuring and executing MCMC simulation 
%   for posterior inference of a Beta Process HMM (BP-HMM) model.
% INPUT:
%  Takes five arguments, each a cell array that specifies parameters
%   as Name/Value pairs, overriding default values (see defaults/ dir)
%  outputParams :
%      specifies how often to save, how often to display, etc.
%      saveEvery, logPrEvery, printEvery, statsEvery, etc.
%  algParams :
%      MCMC alg behavior (# of iterations, parameters of proposal distr.)
%  initParams :
%   initial Markov state (# of initial features, initial state seq., etc.)
%   {'InitFunc', @initBPHMMfromGroundTruth'} initializes to known stateSeq

% Add required libraries, etc.
configNPBayesToolbox;

% ================================================ SANITIZE INPUT
% Converts strings to doubles when possible, etc. to allow command line input
outParams   = sanitizeUserInput( outParams );
algParams   = sanitizeUserInput( algParams );
initParams  = sanitizeUserInput( initParams );

% ================================================ INTERPRET INPUT
algDefs = defaultMCMCParams_BPHMM();
algParams = updateParamsWithUserInput(  algDefs, algParams );

outDefs = defaultOutputParams_BPHMM( outParams, algParams );
outParams = updateParamsWithUserInput( outDefs, outParams(3:end) );

initDefs = defaultInitMCMC_BPHMM();
initParams = updateParamsWithUserInput( initDefs, initParams );

assert( isfield( initParams, 'jobID' ) && isfield( initParams, 'taskID' ) );

INFO = loadSamplerInfo( initParams.jobID, initParams.taskID );
data = INFO.data;
model = INFO.model;
Preproc = INFO.Preproc;


if outParams.doPrintHeaderInfo
    fprintf('Dataset Summary: \n\t %d time series items \n', data.N );
    fprintf('\t   # timesteps per seq:  min %d, median %.1f, max %d ) \n', min(data.Ts), median(data.Ts), max(data.Ts));
    fprintf('\t   %s\n', data.toString() );
end


% ================================================= INITIALIZE MODEL
% Resumes state of previous stored run
%  Psi here is a CHAINHIST struct array, not just a single sampler state
[Psi, algParams, outParams] = initBPHMMResumeRun( data, model, initParams, algParams, outParams );

info_fname = fullfile( outParams.saveDir, 'Info.mat');
save( info_fname, 'data', 'Preproc', 'model', 'initParams', 'algParams', 'outParams' );


% ================================================= RUN INFERENCE
if isfield( algParams, 'TimeLimit' ) && ~isempty( algParams.TimeLimit )
  ChainHist = RunTimedMCMCSimForBPHMM( data, Psi, algParams, outParams, model);
else
  ChainHist = RunMCMCSimForBPHMM( data, Psi, algParams, outParams, model);
end