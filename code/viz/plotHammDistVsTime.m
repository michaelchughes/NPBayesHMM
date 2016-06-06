function [] = plotHammDist( jobIDs, taskIDs, jobNames, objIDs, doHourly )
% View trace plots of hamming distance for sampler state at each
%    recorded iteration of the sampler chain.  Useful for diagnosing
%    convergence and mixing rates.
% ________________________________________________________________________

if ~exist( 'objIDs','var')
    objIDs = 1:6;
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
        DATA = loadSamplerOutput( jobID, taskID , {'times', 'A'} );
        if isnumeric(DATA) && DATA == -1
            continue;
        end
        
        times = DATA.times.Psi;
        
        Ts = DATA.A(end).HDist.Ts( objIDs );
        Hdist = zeros( 1, length(DATA.A ) );
        for ss = 1:length( DATA.A )
            objDists = DATA.A(ss).HDist.obj( objIDs );
            Hdist(ss) = sum( objDists .* Ts )/sum(Ts);
        end
        
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
        
        while length( times ) > 200
           times = times(1:2:end);
           Hdist = Hdist(1:2:end);
        end
        
        if doHourly
            times = times/3600;
        end
        
        plot( times, Hdist, styleStr, 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', taskVis, 'Color', curColor );
        
    end
end

if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    legend( jobNames , 'Location', 'SouthEast' );
else
    legend( taskNames, 'Location', 'SouthEast' );
end

ylabel( 'Hamming Dist', 'FontSize', 18 );

if doHourly
xlabel ('CPU time (hours)', 'FontSize', 18);
else
xlabel ('CPU time (sec)', 'FontSize', 18);
end
grid on;
set( gca, 'FontSize', 16 );
end % main function
