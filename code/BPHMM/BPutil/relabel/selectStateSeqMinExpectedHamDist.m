function [] = selectStateSeqMinExpectedHamDist( jobID, taskIDs, MIN_ITER )
% Take all burned-in samples from multiple chains of MCMC
%   divide into Reference (3/4) and Test (1/4) sets
% Return the stateSeq from the Test set 
%   that minimizes the expected hamming distance to all Reference samples

% jobID = 1121566;
% taskIDs = 1:5;
% MIN_ITER = 15000;
REF_FRAC = 9/10;

% ==================================================== Load Stored Samples
aa = 0;
for tt = 1:length( taskIDs )
    taskID = taskIDs(tt);
    X = loadSamplerOutput( jobID, taskID );
    for ss = 1:length( X.S )
       if X.iters.S(ss) <= MIN_ITER
           continue;
       end
       aa = aa+1;
       AllSamps(aa).stateSeq = X.S(ss).stateSeq;
       AllSamps(aa).taskID = taskID;
       AllSamps(aa).jobID  = jobID;
       AllSamps(aa).sampID = ss;
       AllSamps(aa).sampIter = X.iters.S(ss);
    end
end
    
nSamps = length( AllSamps );
permIDs = randperm( nSamps );
R = round( REF_FRAC * nSamps );
RefSamps = AllSamps( permIDs(1:R) );
TestSamps = AllSamps( permIDs(R+1:end)  );

fprintf( 'Searching for one test sample that best aligns to all reference samples\n');
fprintf( '  |Ref Set|=%d.  |Test Set|=%d\n', R, length(TestSamps)   );

MeanHamDist = zeros( 1, length( TestSamps ) );
tic;
for tt = 1:length( TestSamps )
    HDistToRefs = zeros( 1, R );
    for rr = 1:length( RefSamps )
        HDistToRefs(rr) = calcHammingDistance( RefSamps(rr).stateSeq,  TestSamps(tt).stateSeq );
    end
    MeanHamDist(tt) = mean( HDistToRefs );
    if tt==1 || mod( tt, 25 ) == 0 || tt == length(TestSamps)
        fprintf( ' %5.0f sec. %4d/%d test sequences examined\n', toc, tt, length(TestSamps) );
    end
end

[~, bestID] = min( MeanHamDist );

BestSample = TestSamps( bestID );
X = loadSamplerOutput( BestSample.jobID, BestSample.taskID );
BestSample.A = X.A( BestSample.sampID );

fid = fopen( '~/git/liv-video/SimulationResults.path' );
RESULTS_DIR = textscan( fid, '%s' );
RESULTS_DIR =RESULTS_DIR{1}{1};
fclose(fid);
outpath = fullfile( RESULTS_DIR, num2str(jobID), 'MinExpectedHamDistSample.mat' );
save( outpath, '-struct', 'BestSample' );
fprintf( '... wrote best sample to file %s\n', outpath );