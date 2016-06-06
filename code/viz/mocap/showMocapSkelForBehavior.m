function [] = showMocapSkelForBehavior( jobID, taskID, featIDs, SAVEDIR )
% 
% Uses N. Lawrence's mocap visualization toolbox (download separately)
%  specifically, 
%      acclaimReadSkel.m
%      acclaimLoadChannels.m
%  and skelMultiVis.m

addpath( genpath( '~/MoCap/Erik' ) );

nObjToPlot = 6; % only plot top nObjToPlot objects for each behavior
nSkip = 10;
LEN_fr  =120;
CONTIG_THR   = LEN_fr/12;

[FL, Floc] = getFeatMatHeatmap( jobID, taskID, [], CONTIG_THR );

OUT = loadSamplerOutput( jobID, taskID );
stateSeq = OUT.Psi(end).stateSeq;

INFO = loadSamplerInfo( jobID, taskID );
data = INFO.data;
DATA_DIR = '/data/liv/mhughes/data/MoCapBIG/';

for kk = featIDs
    close all;
    
    [~,rankIDs] = sort( FL(:, kk), 'descend' );
    rankIDs = rankIDs( FL(rankIDs,kk) > 0 );
    
    rankIDs = rankIDs( 1:min( nObjToPlot, length(rankIDs) )  );
    
    for rr = 1:length( rankIDs )    
        seqName = data.seqNames{ rankIDs(rr) };
        
        [~,~,ext] = fileparts( seqName );
        if isempty(ext)
            seqName = [seqName '.amc'];
        end
        
        subjName = [seqName(1:2) '.asf'];
        
        fprintf( ' kk=%d | obj=%d %s\n', kk, rankIDs(rr), seqName );
        ttStart = Floc( rankIDs(rr) ,kk);
        stateSeq( rankIDs(rr) ).z( ttStart:ttStart+9 )
        
        matfilepath = fullfile( DATA_DIR, 'MAT', [seqName '.mat'] );
        if exist( matfilepath, 'file' )
            QQ = load( matfilepath );
            Xch = QQ.Xch;
            skel = QQ.skel;
        else
            skelfilepath = fullfile( DATA_DIR, 'MAT', [subjName '.mat'] );
            if exist( skelfilepath, 'file' )
                QQ = load( matfilepath );
                skel = QQ.skel;
            else
                skel = acclaimReadSkel(  fullfile(DATA_DIR,'asf',subjName) );
            end
            [Xch, skel] = acclaimLoadChannels( fullfile(DATA_DIR,'amc',seqName), skel);
            save( matfilepath, 'Xch', 'skel' );
        end
        
        startT = 12*ttStart;
        curtimes = startT:startT+LEN_fr-1;
        curtimes = curtimes( curtimes <= size(Xch,1) );
        skelMultiVis( skel, Xch( curtimes, : ), 1:nSkip:length(curtimes), [] );

        xlabel(''); ylabel(''); zlabel('');
        set( gca, 'XTick',[],'YTick',[], 'ZTick', []);

        if exist( 'SAVEDIR', 'var' ) && ~isempty( SAVEDIR )
            fname = fullfile( SAVEDIR, sprintf('Behavior%03d-%02d-Seq%s', kk, rr, seqName(1:end-4) ) );
            set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);
           
            export_fig( fname, '-png');
            
            set(gca, 'LooseInset', get(gca, 'TightInset'));

            fnametmp = [fname 'tmp'];
            print(gcf, '-depsc2', '-tiff', '-r300', fnametmp );
            system(['eps2eps ' fnametmp '.eps ' fname '.eps'] );
            system(['rm ' fnametmp '.eps']);

            
%             % Make boundaries tight
%             ti = get(gca,'TightInset')
%             set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% 
%             % Force papersize to be tight as well
%             set(gca,'units','centimeters')
%             pos = get(gca,'Position');
%             ti = get(gca,'TightInset');
% 
%             set(gcf, 'PaperUnits','centimeters');
%             set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             set(gcf, 'PaperPositionMode', 'manual');
%             set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            %export_fig( fname, '-eps'); 
        end
    
    end
    
end