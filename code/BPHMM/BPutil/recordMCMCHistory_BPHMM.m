function ChainHist = recordMCMCHistory_BPHMM( n, outParams, ChainHist, Psi, logPr, Stats )
% Save current state of sampler

% -------------------------------------------------- update logPr trace
if n == 1 || rem( n, outParams.logPrEvery ) == 0
    if isfield( ChainHist, 'logPr' )
        dC = length( ChainHist.logPr ) + 1;
        ChainHist.logPr(dC) = logPr;
    else
        ChainHist = struct();
        ChainHist.logPr = logPr;
        dC = 1;
    end
    
    ChainHist.iters.logPr(dC) = n;
    ChainHist.times.logPr(dC) = toc;
end

% -------------------------------------------------- save current config
if ( n==1 || rem( n, outParams.saveEvery)==0 )
    storePsi = packBPHMMState( Psi );
    if isfield( ChainHist, 'Psi' )
        storeCount = length( ChainHist.Psi ) + 1;
        ChainHist.Psi( storeCount ) = storePsi;
    else
        storeCount = 1;
        ChainHist.Psi = storePsi;
    end
    
    ChainHist.iters.Psi( storeCount) = n;
    ChainHist.times.Psi( storeCount) = toc;
    
    
    ChainHist.RandSeed(storeCount).matlabPRNGState   = RandStream.getGlobalStream.State;
    ChainHist.RandSeed(storeCount).mexPRNGState   = randomseed;
    
end


% -------------------------------------------------- update SamplerAlg stats
if exist( 'Stats','var') && ~isempty( fieldnames( Stats)  )
    if n == 1 || mod( n, outParams.statsEvery ) == 0
        if isfield( ChainHist, 'stats' )
            dC = length( ChainHist.stats ) + 1;
        else
            dC = 1;
        end
        ChainHist.iters.stats(dC) = n;
        ChainHist.times.stats(dC) = toc;
        if isfield( ChainHist, 'TempStats')
            ChainHist = UpdateTempSamplerStats( Stats, ChainHist  );
            ChainHist.stats(dC) = ChainHist.TempStats;
            ChainHist = rmfield( ChainHist, 'TempStats' );
        else
            ChainHist.stats(dC) = Stats;
        end
    else
        ChainHist = UpdateTempSamplerStats( Stats, ChainHist  );
    end
end


end % MAIN FUNCTION


% UpdateTempSamplerStats
%   temporarily stores sampler stats between saves to disk
%    so that we have good records about accept rates throughout the run
function ChainHist = UpdateTempSamplerStats( SamplerStats, ChainHist  )
if ~isfield( ChainHist, 'TempStats' )
    ChainHist.TempStats = SamplerStats;
else
    fNames = fieldnames( SamplerStats );
    for aa = 1:length( fNames )
        Stats = SamplerStats.( fNames{aa} );
        if isfield( Stats, 'C' )
            ChainHist.TempStats.( fNames{aa} ).C = ChainHist.TempStats.( fNames{aa} ).C + Stats.C;
        elseif isfield( Stats, 'nAccept' )
            ChainHist.TempStats.( fNames{aa} ).nAccept = ChainHist.TempStats.( fNames{aa} ).nAccept + Stats.nAccept;
            ChainHist.TempStats.( fNames{aa} ).nTotal = ChainHist.TempStats.( fNames{aa} ).nTotal + Stats.nTotal;
        elseif isfield( Stats, 'ADD' )
            ChainHist.TempStats.( fNames{aa} ).ADD.nAccept = ChainHist.TempStats.( fNames{aa} ).ADD.nAccept + Stats.ADD.nAccept;
            ChainHist.TempStats.( fNames{aa} ).ADD.nTotal = ChainHist.TempStats.( fNames{aa} ).ADD.nTotal + Stats.ADD.nTotal;
            ChainHist.TempStats.( fNames{aa} ).DEL.nAccept = ChainHist.TempStats.( fNames{aa} ).DEL.nAccept + Stats.DEL.nAccept;
            ChainHist.TempStats.( fNames{aa} ).DEL.nTotal = ChainHist.TempStats.( fNames{aa} ).DEL.nTotal + Stats.DEL.nTotal;

        end
    end
end
end
