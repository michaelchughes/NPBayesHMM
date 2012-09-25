function [] = alignStoredSamplerRunToTruth( jobIDs, taskIDs, MIN_ITER )

for jobID = jobIDs
    
    for taskID = taskIDs

        data_struct = loadSamplerInfo( jobID, taskID, '', 'data_struct');
        assert( isfield( data_struct(1), 'true_labels'), 'No ground truth labels to align!' );
            
        X = loadSamplerOutput( jobID, taskID );
        
        meanF = zeros( length( data_struct ), 100 );
        meanFfrac = zeros( length(data_struct), 100);
        nSamps = 0;
        
        tic;
        for storeID = 1:length( X.S )
            
            F = X.S( storeID ).F;
            stateSeq = X.S( storeID ).stateSeq;
            
            % Preprocess F to remove positive entries that 
            %   dont actually explain significant data
            % CUTOFF_THR defines minimize fraction of zseq that 
            %   a positive feature must be assigned in order to count
            Ffrac = double(F);
            for ii = 1:size(F,1)
               ks = find( F(ii,:)==1 );
               for kk = 1:length(ks)
                  kkINDS = stateSeq(ii).z == ks(kk);
                  Ffrac(ii, ks(kk) ) = sum(kkINDS)/data_struct(ii).T;
                  %if ( sum( kkINDS )/data_struct(ii).T ) < CUTOFF_THR
                  %    F(ii, ks(kk) )  = 0;
                  %end
               end
            end
            
            [aF, aFfrac, aStateSeq, Hdist, nTrue] = mapStateSeqLabelsToGroundTruth( F, Ffrac, stateSeq, data_struct );
            
            X.A( storeID ).F = aF;
            X.A( storeID ).Ffrac = aFfrac;
            X.A( storeID ).stateSeq = aStateSeq;
            X.A( storeID ).Hdist = Hdist;
            
            
            if X.iters.S( storeID ) < MIN_ITER
                continue;
            end
            nSamps = nSamps +1;
            meanF( :, 1:size(aF,2) ) = meanF(:,1:size(aF,2) ) + aF;
            meanFfrac(:, 1:size(aF,2) ) = meanFfrac( :, 1:size(aF,2) ) + aFfrac;
        end
        
        lastID = find( sum(meanF,1) > 0, 1, 'last');
        meanF = meanF( :, 1:max(nTrue,lastID) );
        meanF = meanF./nSamps;
        
        meanFfrac = meanFfrac( :, 1:max(nTrue, lastID)  );
        meanFfrac = meanFfrac./nSamps;
        
        X.Summary.meanF = meanF;
        X.Summary.meanFfrac = meanFfrac;
        X.Summary.MIN_ITER  = MIN_ITER;
        
        
        
        fid = fopen( '~/git/liv-video/SimulationResults.path' );
        DATA_DIR = textscan( fid, '%s' );
        fclose(fid);
        DATA_DIR = DATA_DIR{1}{1};
        savefilepath = fullfile( DATA_DIR, num2str(jobID), num2str(taskID), 'SamplerOutput.mat');
        save( savefilepath, '-struct', 'X' );
        
        fprintf( '... %.0f sec | completed align for job %d : %d. saved to file.\n', toc, jobID, taskID );
    end
end


end