function settings = defaultOutputParams_BPHMM( outParams, algP )

saveDir = getUserSpecifiedPath( 'SimulationResults' );

for aa = 1:length( outParams )
    switch aa
        case 1
            jobID = force2double(  outParams{aa} );
        case 2
            taskID = force2double( outParams{aa} );
    end
end
settings.jobID = jobID;
settings.taskID = taskID;
settings.saveDir = fullfile( saveDir, num2str(jobID), num2str(taskID) );
if ~exist( settings.saveDir, 'dir' )
    [~,~] = mkdir( settings.saveDir );
end

if isfield( algP, 'TimeLimit' ) && ~isempty( algP.TimeLimit )
   TL = algP.TimeLimit;
   Niter = Inf;
else
   TL = Inf;
   Niter = algP.Niter;
end

if TL <= 5*60 || Niter <= 200
    settings.saveEvery = 5;
    settings.printEvery = 5;
    settings.logPrEvery = 1;
    settings.statsEvery = 1;
elseif TL <= 2*3600 || Niter <= 5000    
    settings.saveEvery = 25;
    settings.printEvery = 25;
    settings.logPrEvery = 5;    
    settings.statsEvery = 5;
else
    settings.saveEvery = 50;
    settings.printEvery = 50;
    settings.logPrEvery = 10;    
    settings.statsEvery = 10;
end

settings.doPrintHeaderInfo = 1;

%TO DO: Add profiling capability
%settings.Profiler.on           = 0;
%settings.Profiler.baseDir      = getUserSpecifiedPath( 'ProfileResults' );
%settings.Profiler.filename     = 'defaultProfile';
%settings.Profiler.stopAfter    = 5;
