function [] = alignStoredSamplerRunToTruth( jobIDs, taskIDs, BURNFRAC )
% For each given id of a stored sampler runs,
%  align the recovered feature matrix F and state sequences z
%  to the "true" labels available for that data
% For each stored sample in the ChainHistory's Psi field,
%  we add a corresponding "aligned" F,z pair in the 'A' field (A for align)
%  of the ChainHistory
% OUTPUT:
%  no return values
%  Instead, the stored SamplerOutput.mat file is altered
%    to include aligned variables for each stored sampler:
%    A.F     :
%    A.Ffrac :
%    A.stateSeq :

if ~exist('MIN_ITER', 'var')
    BURNFRAC = 0.8;
end

for jobID = jobIDs
    
    for taskID = taskIDs

        data = loadSamplerInfo( jobID, taskID, 'data');
        assert( isprop( data, 'zTrueAll'), 'No ground truth labels to align!' );
            
        ChainHist = loadSamplerOutput( jobID, taskID );
        
        meanF = zeros( data.N, 100 );
        meanFfrac = zeros( data.N, 100);
        nSamps = 0;
        
        tic;
        for storeID = 1:length( ChainHist.Psi )
            
            [alignedPsi] = alignPsiToTruth_OneToOne( ChainHist.Psi(storeID), data );
            
            ChainHist.A( storeID ) = alignedPsi;

            if ChainHist.iters.Psi( storeID ) >= BURNFRAC * ChainHist.iters.Psi( end )
                nSamps = nSamps +1;
                aF = alignedPsi.F;
                meanF( :, 1:size(aF,2) ) = meanF(:,1:size(aF,2) ) + alignedPsi.F;
                meanFfrac(:, 1:size(aF,2) ) = meanFfrac( :, 1:size(aF,2) ) + alignedPsi.Ffrac;
            end

        end
        
        lastID = find( sum(meanF,1) > 0, 1, 'last');
        nTrue = length( unique( data.zTrueAll ) );
        meanF = meanF( :, 1:max(nTrue,lastID) );
        meanF = meanF./nSamps;
        
        meanFfrac = meanFfrac( :, 1:max(nTrue, lastID)  );
        meanFfrac = meanFfrac./nSamps;
        
        ChainHist.Summary.meanF = meanF;
        ChainHist.Summary.meanFfrac = meanFfrac;
        ChainHist.Summary.BURNFRAC  = BURNFRAC;
        
        
        
        fpath = getUserSpecifiedPath( 'SimulationResults' );
        savefilepath = fullfile( fpath, num2str(jobID), num2str(taskID), 'SamplerOutput.mat');
        save( savefilepath, '-struct', 'ChainHist' );
        
        fprintf( '... %.0f sec | completed align for job %d : %d. saved to file.\n', toc, jobID, taskID );
    end
end


end