function [ChainHist, algP, outP] = initBPHMMResumeRun( data, model, initParams, algParams, outParams )

INFO = loadSamplerInfo( initParams.jobID, initParams.taskID );
ChainHist = loadSamplerOutput( initParams.jobID, initParams.taskID );

algP = updateParamsWithUserInput( INFO.algParams, {} );
outP = updateParamsWithUserInput( INFO.outParams, {} );

if isfield( algParams, 'TimeLimit' )
    algP.TimeLimit = algParams.TimeLimit;    
    algP.Niter = [];
else
    algP.Niter = algParams.Niter;
    algP.TimeLimit = [];
end

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


% Remove any entries in the chain history "further in time"
%   than last saved full model state
lastID = find( ChainHist.iters.logPr > ChainHist.iters.Psi(end), 1 );
if ~isempty( lastID )
    ChainHist.iters.logPr = ChainHist.iters.logPr( 1:lastID-1 );
    ChainHist.logPr   = ChainHist.logPr(1:lastID-1);
end
