EXPORT_DIR = '/home/mhughes/git/BPARHMM-NEW/figs/Mocap6/';

Q = load( '/data/liv/mhughes/data/MoCap6/Results/GMM_25trials_K2-20.mat' );
GMM = Q.GMM;
Q = load( '/data/liv/mhughes/data/MoCap6/Results/HMM_25trials_K2-20.mat' );
HMM = Q.HMM;

MNames = {'GMM', 'GMM 1^{st} diff', ...
          'HMM', 'HMM 1^{st} diff', ...
          'BP-HMM'};
MColors = {'b', 'r', 'k'};
MStyles = {'.-', '.--'};
LW = 4;
MS = 30;
figure;
hold all;

Ks = 2:2:20;

for doFirstDiff = [0 1]    
    plot(  Ks,  [ GMM(doFirstDiff+1, Ks).HamDist ] , [MColors{1} MStyles{doFirstDiff+1} ], 'LineWidth', LW, 'MarkerSize', MS  );
end

for doFirstDiff = [0 1]    
    plot(  Ks,  [ HMM(doFirstDiff+1, Ks).HamDist ] , [MColors{2} MStyles{doFirstDiff+1} ] , 'LineWidth', LW, 'MarkerSize', MS  );
end

KKs = 0:2:20;

%jobID = 1121860; taskID = 'best';
jobID = 2134873; taskID = 1;

X = loadSamplerOutput( jobID, taskID );
[~, validIDs, SIDs] = intersect( X.iters.logPr, X.iters.Psi );
[~, bestLid] = max( [ X.logPr(validIDs).all ]  );
bestSid = SIDs(  [X.iters.Psi( SIDs )] == X.iters.logPr( validIDs(bestLid) ) );
BestHamDist = X.A( bestSid ).HDist.total;
plot( KKs, BestHamDist*ones(size(KKs)) ,  [MColors{3}  '--'], 'LineWidth', LW, 'MarkerSize', MS );

% Add Vertical line
plot( 12*ones(1,2), [.15 .75], 'm-.' , 'LineWidth', .75*LW);

set( gca, 'FontSize', 40 );
legend( MNames , 'FontSize', 30);
ylabel( 'Hamming Distance', 'FontSize', 50 );
xlabel( 'Number of Clusters/States', 'FontSize', 50 );
axis( [ 2 20 0.15 .75] );
set( gcf, 'Name', 'HamDistVsK_CompareModels_Mocap6' );
set( gcf, 'Units', 'normalized', 'Position', [0 0 0.5 1] );
set( gcf, 'InvertHardcopy','off','Color',[1 1 1]);

export_fig( fullfile( EXPORT_DIR, ['HamDistVsK_CompareModels_Mocap6']), '-eps');
fprintf( 'Do not forget to manually move the legend so all curves are visible!\n' );