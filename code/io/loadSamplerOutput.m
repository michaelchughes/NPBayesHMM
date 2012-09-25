function [DATA, objIDs] = loadSamplerOutput( jobID, taskID, varNames)
%  Load in save results from sampler run identified by job/task IDs
%  Optional argument varNames is a cell array indicating exactly
%     what fields of the sampler output should be returned.
%  This can be useful if some un-needed fields have very large file size
%     in which case loading everything takes much longer than necessary.

objIDs = [];

DATA_DIR = getUserSpecifiedPath( 'SimulationResults');
DATA_DIR = fullfile( DATA_DIR, num2str(jobID), num2str(taskID) );

dirListing = dir( fullfile( DATA_DIR, 'SamplerOutput.mat' ) );

if isempty( dirListing )
    fprintf( 'Intended sampler output NOT FOUND in %s\n', DATA_DIR );
    DATA = -1;
end

if exist( 'varNames', 'var' )
    needle = '';
    for vc = 1:length( varNames )
       needle = strcat( needle , varNames{vc}, '|');       
    end
    needle = needle(1:end-1); % remove extra pipe delimiter at end    
end

for aa = 1:length( dirListing )
    fullfilename = fullfile( DATA_DIR,  dirListing(aa).name  );
    if exist( 'varNames', 'var' )
        try
            DATA( aa ) = load( fullfilename , '-regexp', needle );
        catch e
            fprintf( 'Corrupt file: %s\n', fullfilename );
            DATA =-1;
            return;
        end
    else
        try
            DATA( aa ) = load( fullfilename );
        catch e
            fprintf( 'Corrupt file: %s\n', fullfilename );
            DATA =-1;
            return;
        end
    end
    
    
    
end


