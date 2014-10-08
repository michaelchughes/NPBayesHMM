clear variables;
EXPORT_DIR = '/home/mhughes/git/BPARHMM-NEW/figs/Mocap6/';

%jobID = 1121860; taskID = 'best';
%jobID = 2134873; taskID = 1;
jobID = 111; taskID = 1;

% ------------------------------------------------------- Load Ground truth
F{1} = zeros(6, 12 );
Q = loadSamplerInfo( jobID, taskID);
for ii = 1:Q.data.N
   activeFeatIDs =  unique( Q.data.zTrue(ii) );
   F{1}(ii, activeFeatIDs) = 1;
end

% ------------------------------------------------------- Load best BPHMM
%  Note that meanFfrac is aggregated across all post-burn-in samples
MIN_THR = 0.02;
% disregard feature if appeared active in less that 2%
%     of state sequence assignment zs when avg'd across all samples
X = loadSamplerOutput( jobID, taskID );
F{2} =  X.Summary.meanFfrac > 0.02;

if size( F{2}, 2 ) > 12
    if all( sum( F{2}(:,13:end), 1 ) == 0 )
        F{2} = F{2}(:,1:12);
    end
end

doTHR = 0;
MINITER = 10000;

keepiters = find( X.iters.Psi >= MINITER );

avgF = zeros( 6, 15 );
for id = keepiters
    K = size( X.A( id ).Ffrac, 2 );
    Fcur = zeros( size(avgF) );
    if doTHR
        Fcur(:,1:K) = X.A( id ).Ffrac > 0.02;
    else
        Fcur(:, 1:K) = X.A( id ).F;
    end
    avgF = avgF + Fcur;
end
avgF = avgF ./ length(keepiters );

kIDs = union( 1:12, find( sum( avgF, 1) > 0) );
F{2} = avgF(:, kIDs);


% ------------------------------------------------------- Load GMM/HMM
Q = load( '/data/liv/mhughes/data/MoCap6/Results/GMM_25trials_K2-20.mat' );
F{3} = Q.GMM(2,12).alignedF;
Q = load( '/data/liv/mhughes/data/MoCap6/Results/HMM_25trials_K2-20.mat' );
F{4} = Q.HMM(2,12).alignedF;

FigNames = {'GroundTruth','BPHMM','GMM','HMM'};
for ff = 1:length( F )
figure;
imagesc( F{ff} );
set( gcf, 'Name', ['Mocap6_FeatMatrix_' FigNames{ff} ]);
colormap( 'bone' );
set( gca, 'FontSize', 20 );

if ff==2 
    rStr = num2str(doTHR);
else
    rStr = '';
end
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);

  export_fig( fullfile( EXPORT_DIR, ['Mocap6_FeatMatrix_' rStr FigNames{ff} ]), '-eps');
end
