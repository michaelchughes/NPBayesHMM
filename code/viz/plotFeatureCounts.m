function [] = plotFeatureCounts( jobIDs, taskIDs,  splitName , jobNames, varName )

if ~exist( 'splitName', 'var' )
    splitName = '';
end
if ~exist( 'varName', 'var' )
    varName = 'total';
end


figure( 'name', 'Feature Counts Visualization'  );
set( gcf, 'Units', 'normalized', 'Position', [0 0.5 0.5 0.5] );
if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    %plotColors  = jet( length( jobNames )  );
    plotColors = get(0,'defaultAxesColorOrder');

    hold on;
else
    plotColors = get(0,'defaultAxesColorOrder');

    hold all;
end

for jobID = jobIDs
    
    for taskID = taskIDs
        DATA = loadSamplerOutput( jobID, taskID, splitName, {'iters', 'S'} );        
        if isnumeric(DATA) && DATA == -1
            continue;
        end
        
        % ---------------------------------  Across All Iterations
        nPerObj = zeros( 3, length( DATA(1).S )  );

        nTotal = zeros( 1, length( DATA(1).S )  );
        nActive = zeros( 1, length( DATA(1).S )  );
        for iter = 1:length( DATA(1).S )
            if strcmp( varName, 'active' )
                nActive(iter) = countActiveStates( DATA.S(iter).F, DATA.S(iter).stateSeq );
            elseif strcmp( varName, 'perObj' )
                nPerObj(1,iter) =  prctile(  sum( DATA.S(iter).F, 2) , 10 );
                nPerObj(2,iter) =  prctile(  sum( DATA.S(iter).F, 2) , 50 );
                nPerObj(3,iter) =  prctile(  sum( DATA.S(iter).F, 2) , 90 );

            end
            
            
            nTotal(iter)  = size( DATA.S(iter).F ,2);
        end
        
        
        if exist( 'jobNames', 'var' )  && ~isempty( jobNames )
            jj = 1 + mod( find( jobID == jobIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            if taskID == taskIDs(1);
                taskVis = 'on';
            else
                taskVis = 'off';
            end
        else
            jj = 1 + mod( find( taskID == taskIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            taskVis = 'on';
        end
        
        if strcmp( varName, 'active' )
            plot( DATA.iters.S,  nActive, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', taskVis, 'Color', curColor );
        elseif strcmp( varName, 'perObj' )
            plot( DATA.iters.S,  nPerObj(2,:) , '-', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', taskVis, 'Color', curColor );
            
            if length( jobIDs ) <= 3
                plot( DATA.iters.S,  nPerObj(1,:) , '--', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', curColor );
                plot( DATA.iters.S,  nPerObj(3,:) , '--', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', curColor );
            end
        else
            plot( DATA.iters.S,  nTotal,  '.-', 'MarkerSize', 15, 'LineWidth', 2, 'HandleVisibility', taskVis, 'Color', curColor );
        end
    end
end

if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    legend( jobNames , 'Location', 'SouthEast' );
else
    taskNames = cell( length(taskIDs), 1 );
    for aa = 1:length( taskIDs )
        taskNames{aa} = num2str( taskIDs(aa) );
    end
    legend( taskNames, 'Location', 'SouthEast' );
end

xlabel ('iteration', 'FontSize', 20);
ylabel( sprintf('%s num. features', varName), 'FontSize', 20 );
set( gca, 'FontSize', 18 );
grid on;
