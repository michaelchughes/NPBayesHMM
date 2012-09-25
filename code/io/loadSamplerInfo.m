function DATA = loadSamplerInfo( jobID, taskID, varNames )
%  Load in original data/settings/model params from sampler run 
%      as identified by job/task IDs
%  Optional argument varNames is a cell array indicating exactly
%     what fields of the sampler output should be returned.
%  This can be useful if some un-needed fields have very large file size
%     in which case loading everything takes much longer than necessary.

DATA_DIR = fullfile( getUserSpecifiedPath('SimulationResults'), ...
                         num2str(jobID), num2str(taskID) );
datafilename = fullfile( DATA_DIR, 'Info.mat' );

if ~exist( 'varNames', 'var' )
    if exist( datafilename, 'file' )
        DATA = load( datafilename );
    else
        error( 'The specified file does not exist.' );
    end
else
    % make sure varNames is a cell arr
    if ischar( varNames )
        varNames = {varNames};
    elseif iscell( varNames )
        % just fine
    else
        error( 'ERROR: expect varNames to be a cell array or a string' );
    end
    
    needle = '';
    for vc = 1:length( varNames )
       needle = strcat( needle , varNames{vc}, '|');       
    end
    needle = needle(1:end-1); % remove extra pipe delimiter at end
     if exist( datafilename, 'file' )
        DATA = load( datafilename, '-regexp', needle );
    else
        error( 'The specified file does not exist.' );
     end
    
    if length( varNames ) == 1
        DATA = DATA.( needle );
    end
     
end

