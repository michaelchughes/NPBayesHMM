function [] = plotSamplerStatsF( jobID, taskID )

DATA = loadSamplerOutput(jobID, taskID);
S = DATA.stats;
iters = DATA.iters.stats;
if isfield( S, 'F' )
if isfield( S(1).F, 'MH' )
    nTry0 = zeros(1,length(S) );
    nAcc0 = zeros(1,length(S) );
    nTry1 = zeros(1,length(S) );
    nAcc1 = zeros(1,length(S) );
    for iterID = 1:length( S );
        C = S(iterID).F.MH.C;
        nTry0(iterID) = sum( C(1,:) );
        nAcc0(iterID) = C(1,2);
        
        nTry1(iterID) = sum( C(2,:) );
        nAcc1(iterID) = C(2,1);
    end
    
    figure;    
    set( gcf, 'Units', 'normalized', 'Position', [0 0 0.5 0.5] );

    set( gcf, 'Name', 'MH Flip Updates of Shared Features' );
    subplot(2, 1, 1);
    plot( iters, nTry1./( nTry1+nTry0 ), 'k--' );
    ylabel( 'Frac. Features On', 'FontSize', 15);
    ylim( [-.01 1.01] );
    grid on;
    set( gca, 'FontSize', 15 );
    
    subplot(2, 1, 2);
    hold on;
    plot( iters, nAcc0./nTry0, 'b.-', 'MarkerSize', 15, 'LineWidth', 2 );
    plot( iters, nAcc1./nTry1, 'r.-', 'MarkerSize', 15, 'LineWidth', 2 );
    legend( '0 --> 1', '1 --> 0' );
    ylabel( 'Frac. Flip Moves Accepted', 'FontSize', 15);
    set( gca, 'YScale', 'log', 'YTick', [1e-9 0.001 0.01 0.1 1], 'YTickLabel', [0 0.001 0.01 0.1 1] )
    set(gca,'XTickMode','manual', 'YTickMode', 'manual');

    grid on; grid(gca,'minor');
    set( gca, 'FontSize', 15 );
end

if isfield( S(1).F, 'RJ' )
    nTryBirth =  zeros(1,length(S) );
    nTryDeath =  zeros(1,length(S) );
    nAccBirth =  zeros(1,length(S) );
    nAccDeath =  zeros(1,length(S) );
    for iterID = 1:length( S );
        RJ = S(iterID).F.RJ;
        nTryBirth(iterID) = RJ.BIRTH.nTrial;
        nAccBirth(iterID) = RJ.BIRTH.nAccept;
        
        nTryDeath(iterID) = RJ.DEATH.nTrial;
        nAccDeath(iterID) = RJ.DEATH.nAccept;
    end
    
    figure;    
    set( gcf, 'Units', 'normalized', 'Position', [0 0.25 0.5 0.5] );
    set( gcf, 'Name', 'RJ Birth-Death Updates of Unique Features' );
    subplot(2, 1, 1);
    plot( iters, nTryBirth./( nTryBirth+nTryDeath ), 'k--' );
    ylabel( 'Frac. Birth Moves', 'FontSize', 15);
    ylim( [-.01 1.01] );
    set( gca, 'FontSize', 15 );
    
    subplot(2, 1, 2);
    hold on;
    birthRate = max(0, nAccBirth./nTryBirth);
    plot( iters, birthRate, 'b.-', 'MarkerSize', 15, 'LineWidth', 2 );
    grid on;

    deathRate = max(0, nAccDeath./nTryDeath);
    plot( iters, deathRate, 'r.-', 'MarkerSize', 15, 'LineWidth', 2 );
    legend( 'Birth', 'Death' );
    ylabel( 'Frac. Moves Accepted', 'FontSize', 15);
    %set( gca, 'YScale', 'log', 'YTick', [1e-9 0.001 0.01 0.1 1], 'YTickLabel', [0 0.001 0.01 0.1 1] )
    ylim( [-.01 1] );
    grid on;
    set( gca, 'FontSize', 15 );
end
end
if isfield( S(1), 'SM' )
    nTryM = zeros(1,length(S) );
    nTryS = zeros(1,length(S) );
    nAccM = zeros(1,length(S) );
    nAccS = zeros(1,length(S) );
    for iterID = 1:length( S );
        SM = S(iterID).SM;
        nTryM(iterID) = SM.MERGE.nTrial;
        nAccM(iterID) = SM.MERGE.nAccept;
        
        nTryS(iterID) = SM.SPLIT.nTrial;
        nAccS(iterID) = SM.SPLIT.nAccept;
    end
    
    figure;
    set( gcf, 'Units', 'normalized', 'Position', [0 0.5 0.5 0.5] );
    set( gcf, 'Name', 'Split-Merge Updates' );
    subplot(2, 1, 1);
    plot( iters, nTryS./( nTryS+nTryM ), 'k--' );
    ylabel( 'Frac. Moves that SPLIT', 'FontSize', 15);
    ylim( [-.01 1.01] );
    set( gca, 'FontSize', 15 );
    
    subplot(2, 1, 2);
    hold on;
    plot( iters, nAccS./nTryS, 'b.-', 'MarkerSize', 15, 'LineWidth', 2 );
    plot( iters, nAccM./nTryM, 'r.-', 'MarkerSize', 15, 'LineWidth', 2 );
    legend( 'Split', 'Merge' );
    ylabel( 'Frac. Moves Accepted', 'FontSize', 15);
    ylim( [0 0.51] );
    %set( gca, 'YScale', 'log', 'YTick', [1e-9 0.001 0.01 0.1 1], 'YTickLabel', [0 0.001 0.01 0.1 1] )
    grid on; 
    set( gca, 'FontSize', 15 );
end

if isfield( S(1), 'U' )
    nTryM = zeros(1,length(S) );
    nTryS = zeros(1,length(S) );
    nAccM = zeros(1,length(S) );
    nAccS = zeros(1,length(S) );
    for iterID = 1:length( S );
        SM = S(iterID).U;
        nTryM(iterID) = SM.BIRTH.nTrial;
        nAccM(iterID) = SM.BIRTH.nAccept;
        
        nTryS(iterID) = SM.DEATH.nTrial;
        nAccS(iterID) = SM.DEATH.nAccept;
    end
    
    figure;
    set( gcf, 'Units', 'normalized', 'Position', [0 0.5 0.5 0.5] );
    set( gcf, 'Name', 'Birth-Death in StateSeq Representation');
    subplot(2, 1, 1);
    plot( iters, nTryS./( nTryS+nTryM ), 'k--' );
    ylabel( 'Frac. Moves that BIRTH', 'FontSize', 15);
    ylim( [-.01 1.01] );
    set( gca, 'FontSize', 15 );
    
    subplot(2, 1, 2);
    hold on;
    plot( iters, nAccS./nTryS, 'b.-', 'MarkerSize', 15, 'LineWidth', 2 );
    plot( iters, nAccM./nTryM, 'r.-', 'MarkerSize', 15, 'LineWidth', 2 );
    legend( 'Death', 'Birth' );
    ylabel( 'Frac. Moves Accepted', 'FontSize', 15);
    ylim( [0 0.51] );
    %set( gca, 'YScale', 'log', 'YTick', [1e-9 0.001 0.01 0.1 1], 'YTickLabel', [0 0.001 0.01 0.1 1] )
    grid on; 
    set( gca, 'FontSize', 15 );
end