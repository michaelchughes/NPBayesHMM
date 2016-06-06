function [F_maxSize, F_maxLoc] = getFeatMatHeatmap( jobID, taskID, objIDs, CONTIG_THR )
% Identify features covering a contig block w/ length at least CONTIG_THR
% Input:
%    jobID : integer id of job
%    taskID : integer id of which run of the job to examine
%    objIDs : vector of objects to consider
%               if empty, we look at all objects
%    CONTIG_THR : integer specifying minimum length of a contig block
%                  mocap data recorded at 120fps,
%                  block-averaged so more like 10 fps
%                  thus setting CONTIG_THR=10 means 1 sec
% OUTPUT
%   F_maxSize : N x K matrix
%                  where entry n,k gives integer length of *longest*
%                   contig block assigned to feature k in sequence n
%   F_maxLoc  : N x K matrix
%                  where entry n,k gives timestep id of the start of the
%                    where *longest* contig block of feat k in seq. n 



if ~exist( 'CONTIG_THR', 'var' )
    CONTIG_THR = 10;
end

INFO = loadSamplerInfo( jobID, taskID );
Q = loadSamplerOutput( jobID, taskID );

data = INFO.data;

F = Q.Psi(end).F;
stateSeq = Q.Psi(end).stateSeq;

K = size(F,2);
if ~exist( 'objIDs', 'var') || isempty( objIDs )
    objIDs = 1:size(F,1);
end
N = length(objIDs);




F_maxSize = zeros( N, K );
F_maxLoc = zeros( N, K );
for aa = 1:length(objIDs )
    [bSiz,bLoc,bVal] = getContigBlocks( stateSeq( objIDs(aa) ).z );
    
    us = unique( bVal );
    for uu = 1:length( us )
        occurIDs = find( bVal == us(uu) );
        
        keepers =  bSiz( occurIDs ) > CONTIG_THR;
        if any( keepers )
            occurIDs = occurIDs(keepers);            
            if isempty( occurIDs )
                continue;
            end
            [maxSize,maxID] = max( bSiz( occurIDs ) );            
            F_maxSize( aa, us(uu) ) = maxSize;
            F_maxLoc( aa, us(uu) ) = bLoc( occurIDs(maxID) );
        end
        
    end
    
end


%           if strcmp( INFO.model.obsModel.type, 'Multinomial' )
%               % Make sure each interval has STIPS in >1/2 of the included tsteps
%               doKeep = true(1, length(occurIDs) );
%               for oo = 1:length( occurIDs )
%                   tstart = bLoc(occurIDs(oo));
%                   tstop  = tstart+bSiz( occurIDs(oo) )-1;
%                   if sum( data_struct(aa).nEmissions(tstart:tstop) ) < length( tstart:tstop )/2
%                       doKeep(oo) = false;
%                   end
%               end
%               occurIDs = occurIDs( doKeep );
%           end