function [] = ReSaveInfoWithGroundTruth( jobID, taskIDs )

data = readSeqDataFromPlainText( '../data/mocap6/' );
data = ARSeqData( 1, data);

for tt = 1:length( taskIDs )
   INFO = loadSamplerInfo( jobID, taskIDs(tt) );
   
   INFO.data = data;
   
   outpath = getUserSpecifiedPath( 'SimulationResults' );
   outpath = fullfile( outpath, num2str(jobID), num2str(taskIDs(tt)), 'Info.mat'  );
   save( outpath, '-struct', 'INFO' );
    
end
