jobIDs = [2130597 2130598 2130599 2130600 2130601 2130602];
taskIDs = 1:4;

infNames = {'SM+zDD AnnealLIN', 'SM+zDD', 'SM+cDD'};
initNames = {'one', 'seq'};

nn=1;
for aa = 1:length(initNames)
for jj = 1:length(infNames)
    jobNames{nn} = [infNames{jj} ' ' initNames{aa}];
    nn=nn+1;
end
end
    
plotLogPrVsTime( jobIDs(1:3), taskIDs, jobNames(1:3) );
set(gcf, 'name', 'OneInit' );
ylim( [-6.45 -6.25]*1e5 );

plotLogPrVsTime( jobIDs(4:6), taskIDs, jobNames(4:6) );
set(gcf, 'name', 'SeqInit' );
ylim( [-6.45 -6.25]*1e5 );
