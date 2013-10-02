% This demo shows examples of data-driven reversible jump (RJ) moves
%  performed on a toy dataset with Gaussian emissions (MANY dims per step)
% We compare and constrast  state-sequence based moves
%   and continuous parameter based moves

clear all;
close all;

doAR = 1;
doDataDriven = 1;
obsDim = 25;

% Generate Gaussian toy data for the BP-HMM
Ktrue = 8;
T = 500;
if doAR
dataP = {'SynthAR', 'nStates', Ktrue, 'obsDim', obsDim, 'T', T};
else
dataP = {'SynthGaussian', 'nStates', Ktrue, 'obsDim', obsDim, 'T', T};
end
mP = {};
initP = {'F.nTotal', 1};
aalgP = {'Niter', 30, ...
         'RJ.birthPropDistr', 'DataDriven', 'BP.doSampleMass', 0, ...
         'BP.doSampleConc', 0, 'HMM.doSampleHypers', 0, ...
          'RJ.minW', obsDim};

for doZMove = [0 1]
    algP = aalgP;
    algP(end+1:end+2) = {'doSampleUniqueZ', doZMove};
    algP(end+1:end+2) = {'doSampleFUnique', ~doZMove};
    
    outP  = {obsDim+(doZMove+1)*100,1};
   
    runBPHMM( dataP, mP, outP, algP, initP );
end

plotLogPr( obsDim+(1:2)*100, 1, {'Unique C', 'Unique Z'} );