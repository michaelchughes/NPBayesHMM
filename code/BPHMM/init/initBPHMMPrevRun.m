function [Psi, algP, outP] = initBPHMMPrevRun( data, model, initParams, algParams, outParams )

INFO = loadSamplerInfo( initParams.jobID, initParams.taskID );
OUT = loadSamplerOutput( initParams.jobID, initParams.taskID );

if isfield( initParams,'queryIter' )
[~, bestID] = min( abs( OUT.iters.Psi - initParams.queryIter ) );
else
    bestID = length( OUT.iters.Psi );
end

Psi = unpackBPHMMState( OUT.Psi( bestID ), INFO.data, INFO.model );

algP = updateParamsWithUserInput( INFO.algParams, {} );
outP = updateParamsWithUserInput( INFO.outParams, {} );
algP.Niter = algParams.Niter;
algP.TimeLimit = algParams.TimeLimit;
outP.saveDir = outParams.saveDir;
outP.jobID   = outParams.jobID;
outP.taskID  = outParams.taskID;

% RESET THE PseudoRandNumGenerator STATE!
mtstream = RandStream('mt19937ar');
RandStream.setGlobalStream(mtstream);

mtstream = RandStream.setGlobalStream(mtstream);
mtstream.State = OUT.RandSeed(bestID).matlabPRNGState;
randomseed(  OUT.RandSeed(bestID).mexPRNGState );
fprintf( 'Initialized from saved job %d %d, set to exact same RNG seed\n', initParams.jobID, initParams.taskID );
