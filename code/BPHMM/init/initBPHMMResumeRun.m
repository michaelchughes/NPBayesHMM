function [ChainHist, algP, outP] = initBPHMMResumeRun( data, model, initParams, algParams, outParams )

INFO = loadSamplerInfo( initParams.jobID, initParams.taskID );
ChainHist = loadSamplerOutput( initParams.jobID, initParams.taskID );

algP = updateParamsWithUserInput( INFO.algParams, {} );
outP = updateParamsWithUserInput( INFO.outParams, {} );
algP.Niter = algParams.Niter;
outP.saveDir = outParams.saveDir;
outP.jobID   = outParams.jobID;
outP.taskID  = outParams.taskID;

% RESET THE PseudoRandNumGenerator STATE!
mtstream = RandStream('mt19937ar');
RandStream.setGlobalStream(mtstream);

mtstream = RandStream.setGlobalStream(mtstream);
mtstream.State = ChainHist.RandSeed(end).matlabPRNGState;
randomseed(  ChainHist.RandSeed(end).mexPRNGState );
fprintf( 'Initialized from saved job %d %d, set to exact same RNG seed\n', initParams.jobID, initParams.taskID );
