% This demo compares several inference algorithms for the BP-HMM
%   on a small toy problem.
% Toy data: Gaussian with  25 sequences, 200 timesteps, and 4 true states
% Inf algs: Prior, Data-Driven, and Split-Merge (SM) moves.
%    see Hughes, Fox, and Sudderth (2012) for algorithm details.
% Each alg is given three random initializations.
%   each of which starts from just one behavior used by all sequences.
% Each chain then runs for 100 iterations.
% We can examine the trace plots of joint log probability
%   to compare the different methods.
%Expected outcome:
% The prior proposal is horribly matched to this data (8 dimensions)
%  so it rarely finds a true state and will likely not recover all 4
% The DD proposal is much better at finding states
%  so it will find all 4 eventually.  Due to a separate proposal for each
%  sequence at each iteration, it will add lots of states in the beginning
% The SM proposal is also very good at finding true states
%   and so it should nicely find all 4 true states (though it is a bit costly)   
% The difference in computation time here highlights
%  why we compared on wall clock time rather than # iterations in the paper

clear variables;
close all;

infAlgs = {'Prior', 'DD', 'SM'};
jobIDs  = [10 11 12];
nTask   = 3;
dataP = {'SynthGaussian', 'nStates', 4, 'nObj', 25, 'T', 200, 'obsDim', 8 };

modelP = {'bpM.gamma', 2}; 
algP_ALL   = {'Niter', 100, 'HMM.doSampleHypers',0,'BP.doSampleMass',0,'BP.doSampleConc',0}; 
initP = {'F.nTotal', 1};
for jj = 1:length(jobIDs)
    jobID = jobIDs(jj);
    
    algP = algP_ALL;
    switch infAlgs{jj}
        case 'Prior'
            algP(end+1:end+2) = {'theta.birthPropDistr', 'Prior'};
        case 'DD'
            algP(end+1:end+2) = {'theta.birthPropDistr', 'DataDriven'};
        case 'SM'
            algP(end+1:end+4) = {'doSplitMerge', 1, 'doSampleFUnique', 0};
    end
   
    for taskID = 1:nTask
        runBPHMM( dataP, modelP, {jobID, taskID, 'printEvery', 25, 'doPrintHeaderInfo', 0}, algP, initP );
    end
    
end

% Plot log prob trace comparison
plotLogPr( jobIDs, 1:nTask, infAlgs );

% Plot recovered theta parameters
figure;
subplot(1,3,1);
plotEmissionParams( jobIDs(1), 1);
title( 'Prior' );
subplot(1,3,2);
plotEmissionParams( jobIDs(2), 1);
title( 'DataDriven' );
subplot(1,3,3);
plotEmissionParams( jobIDs(3), 1);
title( 'Split-Merge' );
