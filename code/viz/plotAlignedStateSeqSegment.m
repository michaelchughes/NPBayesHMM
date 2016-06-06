function [] = plotAlignedStateSeqSegment( jobID, taskID, seqName, timesteps )

INFO = loadSamplerInfo( jobID, taskID );
DATA = loadSamplerOutput( jobID, taskID );

stateSeq = DATA.A(end).stateSeq;

aa = -1;
for ii = 1:length( stateSeq )
    if strcmp( INFO.data.seqNames{ii}, seqName )
        aa = ii;
    end
end
assert( aa > 0, 'Did not find sequence' );
ii = aa;

textSkip = 4;

LetterLabels = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','W','X','Y','Z'};

figure;
imH = imagesc(  stateSeq(ii).z( timesteps ) , [1 26] );
set( imH, 'AlphaData', 0.7 );
set( gca, 'Units', 'normalized', 'Position', [0.1 0.35 0.8 0.6] );
for tt = timesteps(2):textSkip:timesteps(end-1)
        text( tt - timesteps(1)+1, 1, LetterLabels( stateSeq(ii).z(tt) ), 'FontSize', 55, 'HorizontalAlignment', 'Center' );
end
set( gca, 'YTick', [], 'XTick', 1:5:length(timesteps), 'XTickLabel', timesteps(1):5:timesteps(end) );
set( gcf, 'Units', 'Normalized', 'Position', [0 0.7 0.5 0.2] );
set( gca, 'FontSize', 45 );
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);

colormap( [jet(12); spring(14)] );