function [] = plotHammDistPrctile( jobIDs, taskIDs, splitName, jobNames, equivNiters, objIDs )
% View trace plots of hamming distance for sampler state at each
%    recorded iteration of the sampler chain.  Useful for diagnosing
%    convergence and mixing rates.
% ________________________________________________________________________

if ~exist( 'splitName', 'var' )
    splitName = '';
end
if ~exist( 'varName', 'var' )
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

if ~exist( 'objIDs', 'var')
    objIDs = 1:6;
end

for jobID = jobIDs
    taskNames = {};
    
    doFirstTask = 1;
    ttc = 0;
    maxT = 0;
    for taskID = taskIDs
        ttc = ttc+1;
        DATA = loadSamplerOutput( jobID, taskID , splitName, {'iters', 'times', 'S', 'A'}   );
        
        if isnumeric(DATA) && DATA == -1
            continue;
        end
        
        if length( DATA ) == 1
            Ts = DATA.A(end).Hdist.Ts( objIDs );
            
            Hdist = zeros( 1, length(DATA.A ) );
            for ss = 1:length( DATA.A )
                objDists = DATA.A(ss).Hdist.obj( objIDs );
                Hdist(ss) = sum( objDists .* Ts )/sum(Ts);
            end
        else
            error('todo');
        end
        iters = DATA(1).iters.S;
        times = DATA(1).times.S;
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
        
        if ttc == 1
            allHdist = zeros( length(taskIDs), length(times) );
        end

        allHdist(ttc,:) = Hdist; 
        
        curT = find( times > MAX_TIME, 1, 'first' );
        if curT > maxT
            maxT = curT;
        end
        
    end
    
    fprintf( 'Job %d using iters 0-%d\n', jobID, DATA.iters.S(maxT ));
    
    HD50 = prctile( allHdist ,50 );
    HD05 = prctile( allHdist , 10 );    
    HD95 = prctile( allHdist , 90 );
    
    ts = linspace( 0, MAX_TIME, maxT );
    
    plot( ts, HD50(1:maxT), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', 'on', 'Color', curColor );
    plot( ts, HD05(1:maxT), '--', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', curColor );
    plot( ts, HD95(1:maxT), '--', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', curColor );
end

if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    legend( jobNames , 'Location', 'SouthEast' );
else
    legend( taskNames, 'Location', 'SouthEast' );
end

if strcmp( varName, 'all' )
    varName = '';
end

ylabel( 'Hamming Dist', 'FontSize', 18 );


if isfield( DATA,'times')
    xlabel ('cpu time (sec)', 'FontSize', 18);
else
    xlabel ('iteration', 'FontSize', 18);
end

grid on;
set( gca, 'FontSize', 16 );
end % main function
