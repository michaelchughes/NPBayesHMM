function [] = plotStateSeq( varargin )
% Visualize hidden state seq of the HMM fit to given collection of
%   sequential data. Will plot Ground Truth labels as well if available.
%USAGE:
%   Loading from MAT file, plotting specific objects by ID
%         plotStateSeq( jobID, taskID,  queryIter, objIDs  )
%   Showing the state sequence already loaded in memory as workspace var
%         plotStateSeq( stateSeq, objIDs )
%CAVEATS:
%   Can only display 12 sequences at once (otherwise looks too small)
%    also, will only display first 5000 timesteps of each sequences
%   To change these settings, use parameters defined right below.

MAX_NUM_OBJ = 12;
Tmax = 5000;

GT = [];
% ================================================ Process User Input
if isstruct( varargin{1} )
    if isfield( varargin{1}, 'stateSeq' )
        stateSeq = varargin{1}.stateSeq;
    elseif isfield( varargin{1}, 'z' )
        stateSeq = varargin{1};
    end
    
    for ll = 1:length(varargin)
       if isobject( varargin{ll} )
         GT = varargin{ll};
       elseif isnumeric( varargin{ll} )
         objIDs = varargin{ll};
       end
    end
else
    jobID = varargin{1};
    taskID = varargin{2};
    DATA = loadSamplerOutput( jobID, taskID, {'iters', 'Psi'} );
    data = loadSamplerInfo( jobID, taskID, {'data'} );

    if isprop( data, 'zTrueAll') && ~isempty( data.zTrueAll )
        GT = data;
    end
    
    if length( varargin ) >= 3 && varargin{3} >= 0
        queryIter = varargin{3};
        [~, idx] = min( abs( queryIter - DATA.iters.Psi ) );
    else
        idx = length( DATA.Psi );
    end
    Psi = DATA.Psi( idx );
    stateSeq = Psi.stateSeq;
    
    if length( varargin ) >= 4
        objIDs = varargin{4};
    end
end

if ~exist( 'objIDs', 'var' )
    objIDs = 1:length( stateSeq );
end

if length( objIDs ) > MAX_NUM_OBJ
    objIDs = objIDs( 1:MAX_NUM_OBJ );
end
fprintf( 'Showing sequences %d - %d\n', objIDs(1), objIDs(end) );
nObj = min( MAX_NUM_OBJ, length( objIDs ) );

% ------------------------------------------  Read in first Tmax of each sequence
Ts = zeros( nObj, 1 );
for aa = 1:nObj
    Ts(aa) = min( Tmax, length( stateSeq( objIDs(aa) ).z )   );
end
aggT = [0; cumsum( Ts )];


zEst = zeros( 1, aggT(end) );
if ~isempty(GT)
    zTrue = zeros(1, aggT(end) );
end
for aa=1:nObj
    if isfield( stateSeq(1), 'zHat' )
        zEst( aggT(aa)+1:aggT(aa+1)) = stateSeq(  objIDs(aa)    ).zHat(1:Ts(aa) );
    else
        zEst( aggT(aa)+1:aggT(aa+1)) = stateSeq(  objIDs(aa)    ).z(1:Ts(aa) );
    end
    if ~isempty( GT )
        zzTrue_aa = GT.zTrue( objIDs(aa) );
        zTrue( aggT(aa)+1:aggT(aa+1) ) = zzTrue_aa( 1:Ts(aa) );
    end
end

if exist( 'zTrue', 'var' )
    [zEstR, TrueAlphaNames, EstAlphaNames] = alignEstStateSeqToTruth( zEst, zTrue );
    uEstR = unique(zEstR);
    uTrue = unique( zTrue );
    MyColors = jet( 5*length( uTrue ) );
    LabelRange = [min(unique(zEstR)) max(uTrue)];
else
    zEstR = zEst;
    uEstR = unique(zEstR);
    MyColors = jet( length( uEstR ) );
    LabelRange = [min(uEstR) max(uEstR)];
end

if LabelRange(1)==LabelRange(2)
   LabelRange(1)=LabelRange(1)-.01; 
end


% ============================================================ PLOT
halfV = 5; % defines how "tall" each seq's plot is... large allows real estate for interaction
aggii = 0;
hh = figure( );

set(hh,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
set( hh,'Visible', 'Off' );
drawnow;

nRows = length( objIDs ); % num rows for sub plot
nCols = 1;

aa = 0;
for ii = objIDs
    aa = aa + 1;
    aggii = aggii + 1;
    T_ii    = Ts( aa );
    
    subplot(nRows,nCols,aa);
    hold on;
    
    origzImg = zEst( aggT(aggii)+1:aggT( aggii+1 )  );
    curzImg = zEstR( aggT(aggii)+1:aggT( aggii+1 )  );
    if ~isempty( GT )
        curzTrue = zTrue( aggT(aggii)+1:aggT( aggii+1 )  );
        
        % Plot recovered and true seq. one on top of the other
        stateSeqImg_ii = [ repmat( curzImg,    halfV, 1); ...
                           repmat( curzTrue, halfV, 1) ];
        
        h = imagesc(stateSeqImg_ii, LabelRange);
        
        colormap( MyColors );
        
        % Draw thin line to separate
        plot( 1:T_ii, 0.5+halfV*ones(T_ii,1), 'k--', 'LineWidth', 0.5 );
        
        % Annotate with TRUE state labels
        [bS, bPos, bLabel] = getContigBlocks( curzTrue );
        for bb = 1:length( bLabel )
            if bS( bb ) > 10
                alphaID = uTrue == bLabel(bb);
                text( bPos(bb)+0.5*bS(bb), 1.5*halfV, TrueAlphaNames{alphaID}, 'FontSize', 12 , 'HorizontalAlignment','center');
            end
        end
        
        % Annotate with EST state labels
        uEst = unique( zEst );
        [bS, bPos, bLabel] = getContigBlocks( origzImg );
        for bb = 1:length( bLabel )
            if bS( bb ) > 10
                alphaID = uEst == bLabel(bb);
                text( bPos(bb)+0.5*bS(bb), 0.5*halfV, EstAlphaNames{alphaID}, 'FontSize', 12, 'HorizontalAlignment','center' );
            end
        end
        
        % Annotate TRUE and ESTIMATED rows
        text( -15, 1.5*halfV, 'TRUE' );
        text( -15, 0.5*halfV, 'EST' );
        
    else
        stateSeqImg_ii = repmat( curzImg, 2*halfV, 1 );
        h = imagesc(stateSeqImg_ii, LabelRange );
        colormap( MyColors );
    end
    
    set(  h, 'AlphaData', 0.5 );
    
    % Make plot look nice:
    %      remove y labels (unneccessary for interpreting)
    %      adjust axes to remove whitespace
    %      adjust length scales to be roughly similar
    axPos = get( gca, 'Position' );
    %maxT = prctile( Ts, 80 );
    maxT = max( Ts );
    set( gca, 'Position', [axPos(1) axPos(2)  axPos(3)*min(1,T_ii/maxT )  0.8*axPos(4) ] )
    set( gca, 'YTick', [] );
    
    set( gca,'CLim', LabelRange );
    xlim(gca,  [ 0.5 Ts(aggii)+0.5 ] );
    ylim(gca,  [ 1 2*halfV] );
    
    
end % loop over time series in current category

set( hh,'Visible', 'On' );

drawnow;
end % main function



