close all;
jobNames = {'SM+cDD', 'SM+cDD AnnealEXP', 'SM+zDD', 'SM+zDD AnnealLIN', 'SM+cDD no H'};
jobIDs = [1515157 1876679 1878190 1878692 1515158];

plotLogPrVsTime( jobIDs(1:4), 1:10, jobNames(1:4) );
ylim( [-5.55 -5.35]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_cDD&zDD' );

plotLogPrVsTime( jobIDs(1:2), 1:10, jobNames(1:2) );
ylim( [-5.55 -5.35]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_cDDOnly' );

plotLogPrVsTime( jobIDs([1:2 5]), 1:10, jobNames([1:2 5]) );
ylim( [-5.55 -5.35]*1e4 ); xlim( [0 35000]);
set(gcf,'Name', 'LogPr_cDDOnly++' );

for jj = 1:length( jobNames )
plotStateSeq(  jobIDs(jj), 'best' );
set(gcf, 'Name', sprintf('StateSeq_best_%s', jobNames{jj})  );
end

