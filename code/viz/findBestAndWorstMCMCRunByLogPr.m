function [] = findBestAndWorstMCMCRunByLogPr( jobIDs, nTasks, BURN_FRAC )
% Given list of jobIDs (each with many runs of MCMC)
%   identify for each one the BEST and WORST runs by log probability
% Save these as symbolic directories
%      <path/to/results>/<jobID>/best/
%      <path/to/results>/<jobID>/worst/
% INPUT
%   jobIDs : row vector of job IDs
%   nTasks : row vector where entry jj says # of tasks in job jobIDs(jj)
%               if scalar, assume same # tasks for all job IDs
%   BURN_FRAC : scalar fraction of samples to discard when analyzing 
%                log prob. 
%           Example: BURN_FRAC=4/5 means only retain last 1/5 of samples

if ~exist( 'BURN_FRAC', 'var' )
    BURN_FRAC = 4/5;
end

if length( nTasks )==1 && length(jobIDs) > 1
   nTasks = repmat( nTasks, 1, length(jobIDs) ); 
end

RESULTS_DIR = getUserSpecifiedPath( 'SimulationResults' );

for jj = 1:length( jobIDs ); 
    
   jobID = jobIDs(jj);
    
   meanLogPr = NaN(1, nTasks(jj) );
   for taskID = 1:nTasks(jj)
       OUT = loadSamplerOutput( jobID, taskID , {'iters', 'logPr'}   );
       if ~isnumeric( OUT )
           lastpt = length( OUT.logPr );
           startpt = ceil( BURN_FRAC*lastpt);
           meanLogPr( taskID ) = mean( [OUT.logPr( startpt:lastpt ).all ] );
       end
   end
        
   [~, bestTaskID] = max( meanLogPr );
   [~, worstTaskID] = min( meanLogPr );
   
   curDir = pwd;
   cd( fullfile( RESULTS_DIR, num2str( jobID ) )  );
   
   targetDir = [num2str( bestTaskID ) '/' ];
   linkDir = ['./best'];
   if exist( linkDir, 'file' )
       system( ['rm ' linkDir ] );
   end
   MAKE_LINK_CMD = ['ln -s ' targetDir ' ' linkDir ];
   [status, result] = system( MAKE_LINK_CMD ); 
   
   targetDir = [num2str( worstTaskID ) '/' ];
   linkDir = ['./worst'];
   if exist( linkDir, 'file' )
       system( ['rm ' linkDir ] );
   end
   MAKE_LINK_CMD = ['ln -s ' targetDir ' ' linkDir ];
   [status, result] = system( MAKE_LINK_CMD ); 
   
   cd( curDir );
end