addpath(genpath( '~/code/exportfig/') );

close all;
jobNames = {'SM+zDD+Anneal', 'SM+zDD', 'SM+cDD+Anneal', 'SM+cDD', 'Prior Rev. Jump'};
jobIDs = [2134873 2134874 2134875 2134876 2207266] %1515159];

fSize = 30;
objIDs = 1:6;
doHourly = 1;

EXPORT_DIR = '/home/mhughes/git/BPARHMM-NEW/figs/Mocap6/';

plotLogPrVsTime( jobIDs([1 2 4 5]), 1:25, jobNames([1 2 4 5]), doHourly );
ylim( [-5.8 -5.34]*1e4 );
xlim( [0 10]);
grid off;
set( gcf, 'units', 'normalized', 'position', [0 0.5, 0.5 0.5] );
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);
set(gcf,'Name', 'Mocap6_LogPr_CompareAnnealingInference' );
set( gca, 'FontSize', fSize );
xlabel( 'CPU time (hours)', 'FontSize', fSize+3 );
ylabel( 'joint log prob.', 'FontSize', fSize+3 );
legend('off');
%[legh,objh,outh,outm] = legend;
%set( objh, 'LineWidth', 4 );
%set( legh, 'FontSize', fSize+3 );
%export_fig( fullfile(EXPORT_DIR, 'Mocap6_CompareInfer_LogPr'), '-eps');

plotHammDistVsTime( jobIDs([1 2 4 5]), 1:25, jobNames([1 2 4 5]), objIDs, doHourly );
ylim( [0.1 0.6] );
xlim( [0 10]);
grid off;
set( gcf, 'units', 'normalized', 'position', [0 0.5, 0.5 0.5] );
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);
set(gcf,'Name', 'Mocap6_HamDist_CompareAnnealingInference' );
legend('off');
%legend( 'Location', 'EastOutside');
set( gca, 'FontSize', fSize );
xlabel( 'CPU time (hours)', 'FontSize', fSize+3 );
ylabel( 'Hamming dist.', 'FontSize', fSize+3 );
%[legh,objh,outh,outm] = legend;
%set( objh, 'LineWidth', 4 );
%set( legh, 'FontSize', fSize+3 );
%export_fig( fullfile(EXPORT_DIR, 'Mocap6_CompareInfer_HamDist'), '-eps');


% Make common legend
plotHammDistVsTime( jobIDs([1 2 4 5]), 1, jobNames([1 2 4 5]), [], doHourly );
set( gcf, 'units', 'normalized', 'position', [0 .5 0.2 0.2] );
legend( 'Location', 'East');
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);
[legh,objh,outh,outm] = legend;
set( objh, 'LineWidth', 6 );
set( legh, 'FontSize', fSize+5 );
set( gca, 'Visible', 'off' );
%export_fig( fullfile(EXPORT_DIR, 'Mocap6_CompareInfer_legend'), '-eps');