close all;
jobNames = {'SM+zDD+AnnealLin', 'SM+zDD', 'SM+cDD+AnnealLin', 'SM+cDD'};
jobIDs = [2134873 2134874 2134875 2134876 ];

for jj = jobIDs
    findBestAndWorstMCMCRunByLogPr(jj, 25);
end

plotLogPrVsTime( jobIDs(1:2), 1:25, jobNames(1:2) );
ylim( [-5.55 -5.34]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_zDD' );


plotLogPrVsTime( jobIDs(3:4), 1:25, jobNames(3:4) );
ylim( [-5.55 -5.34]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_cDD' );


plotLogPrVsTime( jobIDs(1:4), 1:25, jobNames(1:4) );
ylim( [-5.55 -5.34]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_CompareAnnealingInference' );


for jj = 1:length( jobNames )
    plotStateSeq(  jobIDs(jj), 'best' );
    set(gcf, 'Name', sprintf('StateSeq_best_%s', jobNames{jj})  );
end
