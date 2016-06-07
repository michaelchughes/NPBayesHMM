function [] = plotFeatMat( varargin )
% Show the binary feature matrix of current sampler state.
% Uses "bone" colormap, so *active* features (F==1) are WHITE, F==0 black.
% Designed to handle diverse inputs easily. When an HMM state sequence 
%   is provided, also shows difference between *active* and *available*,
%   where both imply F(ii,kk)==1, but *active* additionally requires that
%   at least one sequence has at least one timestep in z assigned to kk.
%USAGE:  
%  To see matrix for current matrix "F"
%    plotFeatMat( F )
%  To see matrix for current whole model config "Psi" (where F=Psi.F)
%    plotFeatMat( Psi )
%  To see matrix for stored sampler run, at particular iteration
%    plotFeatMat( jobID, taskID, queryIter )

if isstruct( varargin{1} )
    if isfield( varargin{1}, 'F' )
        F = varargin{1}.F;
        stateSeq = varargin{1}.stateSeq;
    end
    
    for ll = 1:length(varargin)
       if isnumeric( varargin{ll} )
         objIDs = varargin{ll};
       end
    end
elseif (isnumeric( varargin{1} ) || islogical(varargin{1}) )&& size( varargin{1}, 1) > 1
    F = varargin{1};
else
    jobID = varargin{1};
    taskID = varargin{2};
    DATA = loadSamplerOutput( jobID, taskID, {'iters', 'Psi'} );
       
    if length( varargin ) >= 3 && varargin{3} >= 0
        queryIter = varargin{3};
        [~, idx] = min( abs( queryIter - DATA.iters.Psi ) );
        fprintf( '@ iter %d\n', DATA.iters.Psi(idx) );
    else
        idx = length( DATA.Psi );
    end
    Psi = DATA.Psi( idx );
    F = Psi.F;
    stateSeq = Psi.stateSeq;
    
    if length( varargin ) >= 4
        objIDs = varargin{4};
    end
end


if exist('stateSeq', 'var')
    F_used = zeros(size(F));
for ii = 1:size(F,1)
    F_used(ii,unique( stateSeq(ii).z)) = 1;
end
else
    F_used = F;
end

if ~exist( 'objIDs', 'var' )
    objIDs = 1:size(F,1);
else
    fprintf( 'Showing sequences %d - %d\n', objIDs(1), objIDs(end) );
end
F = F(objIDs,:);
F_used = F_used(objIDs,:);

imagesc(F+F_used, [0 2]);

xlabel( 'behaviors', 'FontSize', 20);
ylabel( 'data sequences', 'FontSize', 20 );
set( gca, 'FontSize', 16 );
colormap bone;
colorbar(  'YTick', [0 1 2],       'YTickLabel', {'Disabled', 'Available', 'Active'});

% Force ytick labels to be integers
ticks = get( gca, 'YTick');
ticks = ticks( ticks >= 1 );
ticks = ticks( ticks <= length(objIDs) );
ticks = unique(int32(ticks));
set( gca, 'YTickLabel', objIDs(ticks) );

end