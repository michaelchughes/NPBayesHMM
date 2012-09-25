function settings = defaultOutputParams_BPHMM( outParams, Niter )

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
if Niter <= 200
    settings.saveEvery = 5;
    settings.printEvery = 5;
    settings.logPrEvery = 1;
    settings.statsEvery = 1;
elseif Niter <= 5000    
    settings.saveEvery = 25;
    settings.printEvery = 50;
    settings.logPrEvery = 5;    
    settings.statsEvery = 5;
else
    settings.saveEvery = 50;
    settings.printEvery = 100;
    settings.logPrEvery = 10;    
    settings.statsEvery = 10;
end

settings.Profiler.on           = 0;
settings.Profiler.baseDir      = getUserSpecifiedPath( 'ProfileResults' );
settings.Profiler.filename     = 'defaultProfile';
settings.Profiler.stopAfter    = 5;

settings.doPrintHeaderInfo = 1;
