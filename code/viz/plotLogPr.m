function [] = plotLogPr( jobIDs, taskIDs, jobNames, varName )
% View trace plots of joint log probability for sampler state at each
%    recorded iteration of the sampler chain. Useful for diagnosing
%    convergence and mixing rates, and comparing multiple runs.
% SYNTAX:  plotLogPr( jobIDs, taskIDs, jobNames*, varName* ) where *=optional
% USAGE:
%    To view a single sampler run's log probability trace, 
%       plotLogPr( <jobID>, <taskID> )
%    To view log prob. of samples for a particular variable named "X"
%       plotLogPr( <jobID>, <taskID>, {}, 'X' )
%    To compare results from multiple jobs
%       plotLogPr( [jobA, jobB], taskIDs, {'A', 'B'} )
%    To compare multiple runs for multiple jobs
%       plotLogPr( <jobID vector>, <taskID vector>, {'A', 'B', 'C', ...}  )
% ______________________________________________________________________________

if ~exist( 'varName', 'var' ) || isempty( varName )
    varName = 'all';
end

figure;
set( gcf, 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5] );
if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    plotColors  = get(0,'defaultAxesColorOrder');
    hold on;
else
    plotColors = get(0,'defaultAxesColorOrder');
    hold all;
end


for jobID = jobIDs
    taskNames = {};
    
    doFirstTask = 1;
    for taskID = taskIDs
        DATA = loadSamplerOutput( jobID, taskID , {'iters', 'times', 'logPr'} );
        if isnumeric(DATA) && DATA == -1
            continue;
        end
        
        logPr = [ DATA.logPr(:).( varName ) ];
       
        iters = DATA.iters.logPr;

        if exist( 'jobNames', 'var' )  && ~isempty( jobNames )
            jj = 1 + mod( find( jobID == jobIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            if doFirstTask;
                taskVis = 'on';
                doFirstTask = 0;
            else
                taskVis = 'off';
            end
        else
            taskNames{end+1} = num2str( taskID );

            jj = 1 + mod( find( taskID == taskIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            taskVis = 'on';
        end
        
        styleStr = '.-';
        plot( iters, logPr, styleStr, 'MarkerSize', 15, 'LineWidth', 2, ...
               'HandleVisibility', taskVis, 'Color', curColor );
    end
end

if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    legend( jobNames , 'Location', 'SouthEast' );
else
    legend( taskNames, 'Location', 'SouthEast' );
end

if strcmp( varName, 'all' )
    varName = '';
end

ylabel( sprintf('log prob. %s', varName), 'FontSize', 18 );

xlabel ('iteration', 'FontSize', 18);

grid on;
set( gca, 'FontSize', 16 );
end % main function
