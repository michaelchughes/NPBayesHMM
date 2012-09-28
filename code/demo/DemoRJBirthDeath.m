% This demo shows examples of data-driven reversible jump (RJ) moves
%  performed on a toy dataset with Gaussian emissions (6 dims per timestep)
% To focus on the proposal procedure, we only use ONE sequence
%  and show several iterations of the RJ moves, interleaved with
%  other Gibbs updates for the already existing emission parameters.
%  Of course, this demo is just a subset of all possible sampler moves.
% At each iteration, we show the proposed theta parameters in *red*
%   alongside existing params (blues). We also show the data sequence
%   and (if DD proposals are used) show the window that yielded the
%   proposals, colored according to accept/rejection of the proposal.
% Under DD moves, the proposal distribution is a *mixture* of the prior
%  and the posterior in the subwindow, so we have 1/2 chance of drawing
%  a prior proposal, and 1/2 chance drawing from the window.
%  We mark clearly which was used for each step as either "prior" or "DD".
% Note that the actual observations have *six* dimensions,
%  and the inference is performed on all six dimensions,
%  but for visualization we restrict ourselves to the first 2
% To examine behavior with the DD proposals vs. just the Prior propsals,
%  simply toggle the "doDataDriven" flag from 1 (=DD) to 0 (=Prior)

clear all;
close all;

doDataDriven = 1;

% Generate Gaussian toy data for the BP-HMM
%   4 true behaviors, each observation is a 20-dimensional vector
%   only 1 sequence of 1000 timesteps (T=1000)
Ktrue = 4;
if ~exist('data','var')
    [data, True] = genToySeqData_Gaussian( Ktrue, 6, 1, 1000, 1 );
end

mP = defaultModelParams_BPHMM(data);
mP.bpM.gamma = 2;

initP = defaultInitMCMC_BPHMM();
algP  = defaultMCMCParams_BPHMM();
outP  = defaultOutputParams_BPHMM({1,1}, algP);

if doDataDriven
    algP.theta.birthPropDistr = 'DataDriven';
else
    algP.theta.birthPropDistr = 'Prior';
end

initP.F.nTotal = 1;
outP.doPrintHeaderInfo = 0;
Psi = initBPHMMFresh( data, mP, initP, algP, outP);

figure(1); clf;
set( gcf, 'units', 'normalized', 'position', [.25 0.6 0.8 0.3] );
plotData( data, 1 );

figure(101); clf;
set(gcf, 'units', 'normalized', 'position', [.6 .1 .3 .5] );
hold all;
plotEmissionParams(Psi, data);
title( 'theta (@ init)' );
drawnow;
pause;

for trial = 1:40
    
   [Psi, Stats, R] = sampleUniqueFeats(Psi, data, algP );
   
   figure(1); clf;
   plotData(data,1);
   if R.doAccept
       resultStr = 'ACCEPT';
       colVec = [0 1 0.5];
   else
       resultStr = 'reject';
       colVec = [0.7 0 0.1];
   end
   if R.doBirth       
       if R.choice == 1
            fprintf('%s | birth feat %d | prior \n', resultStr, R.activeFeatIDs );
       else
            fprintf('%s | birth feat %d | DD \n', resultStr, R.activeFeatIDs );

            xs = R.window(1):R.window(2);
            plot( xs, 2*ones(size(xs)), '-', 'Color', colVec, 'LineWidth', 5 );
            plot( xs, -2*ones(size(xs)), '-', 'Color', colVec, 'LineWidth', 5); 
            plot( xs(1)*[1 1], [-2 2], '-', 'Color', colVec, 'LineWidth', 5);
            plot( xs(end)*[1 1], [-2 2], '-', 'Color', colVec, 'LineWidth', 5);
       end
       
   else
      fprintf('%s | death feat %d\n', resultStr, R.activeFeatIDs ); 
   end
   
   figure(101); 
   plotGauss( R.thetaStar.mu', R.thetaStar.invSigma, 1, [1 0 0] );
   drawnow;
   pause(0.5);
   
   Psi = sampleStateSeq( Psi, data );
   Psi.ThetaM = Psi.ThetaM.sampleAllTheta( data, Psi.stateSeq );
   Psi.TransM = Psi.TransM.sampleAllEta( Psi.F, Psi.stateSeq );
   
   figure(101); clf;
   hold all;
   plotEmissionParams(Psi, data);
   title( ['theta @iter' num2str(trial)] , 'FontSize', 20);
   drawnow;
   
end