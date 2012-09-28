close all;
%findBestAndWorstMCMCRunByLogPr( 1515157:1515161, 10 );

% plotLogPrVsTime( [1515157 1515159 1515160], 1:10, {'SM+DD one', 'Prior u5', 'SM+DD u5'} )
% ylim( [-5.8 -5.33]*1e4 );
% set( gcf, 'Name', 'LogPrCompare_ValidSamplers' );
% set(gca,'FontSize', 20);
% 
% plotLogPrVsTime( [1515157 1515158], 1:10, {'SM+DD one', 'SMnoq+DD one'} )
% ylim( [-5.8 -5.33]*1e4 );
% set( gcf, 'Name', 'LogPrCompare_HastingsFactor_ONE' );
% set(gca,'FontSize', 20);
% 
% plotLogPrVsTime( [1515160 1515161], 1:10, {'SM+DD unique5', 'SMnoq+DD unique5'} )
% ylim( [-5.8 -5.33]*1e4 );
% set( gcf, 'Name', 'LogPrCompare_HastingsFactor_Unique5' );
% set(gca,'FontSize', 20);


plotStateSeq( 1515157, 'best' );
set(gcf,'Name', 'StateSeq_SMvalid_ONE_Best' );

plotStateSeq( 1515158, 'best' );
set(gcf,'Name', 'StateSeq_SMnoq_ONE_Best' );

plotStateSeq( 1515160, 'best' );
set(gcf,'Name', 'StateSeq_SMvalid_UNIQUE5_Best' );

plotStateSeq( 1515161, 'best' );
set(gcf,'Name', 'StateSeq_SMnoq_UNIQUE5_Best' );
