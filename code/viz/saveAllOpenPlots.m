function [] = saveAllOpenPlots( dirName, figExt )

MY_DIR = getUserSpecifiedPath( 'Figures' );

if ~exist('figExt', 'var' )
    figExt = 'png';
end

figHandles = get( 0, 'Children' );

for hh = figHandles'
    figName = get( hh, 'Name' );
        
    if isempty( figName)
        continue;
    end
    
    [~,figName,~] = fileparts(figName);
    savefilename = sprintf( '%s.%s', figName, figExt  );
    
    saveDir = fullfile( MY_DIR, dirName);
    [~,~] = mkdir( saveDir );
    
    savefilename = fullfile( saveDir, savefilename );
    switch figExt
        case 'eps'
            figure(hh);
            set(hh, 'PaperPositionMode', 'auto');
            print(savefilename,'-depsc2','-r300');
            
        case 'png'
            set(hh, 'PaperPositionMode', 'auto');
            saveas( hh, savefilename, figExt );
            
        case 'fig'
            %set(hh, 'PaperPositionMode', 'auto');
            saveas( hh, savefilename, figExt );
    end
    fprintf( '.....................  saved to file %s\n', savefilename );
end