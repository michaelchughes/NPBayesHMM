% readJointAnglesAsMatrixFromAMC.m
%  Read data from AMC motion capture file into Matlab matrix
% CREDITS
%  Modified from amc_to_matrix.m (Jernej Barbic, CMU, March 2003)
%    with further changes by E.Sudderth and E. Fox (MIT, 2008-2009)
%  Smoothing Angle suggestion due to N. Lawrence (smoothAngleChannels.m)

function [D,ColNames] = readJointAnglesAsMatrixFromAMC( fname, QueryChannelNames )

% Preallocate the matrix
%  since we usually have HUGE volumes of sensor data
D = nan( 7200, 50 );

% Open file
fid=fopen(fname, 'rt');
if (fid == -1)
  error('ERROR: Cannot open file %s.\n', fname);
end;

% Read lines until we've skipped past the header 
%   assuming it ends with the line ":DEGREES"
line=fgetl(fid);
while ~strcmp(line,':DEGREES')
  line=fgetl(fid);
end

% Loop through each frame, one-at-a-time
fID = 0;
%partID = 0;
%nDimsPerPart = [];
fNames = {};
fCtr = 1;
while ~feof(fid)
    
  line = fgets( fid );
  SpaceLocs = strfind( line, ' ');
  
  if isempty( SpaceLocs ) && ~isempty(line)
      % Advance to next time frame (fID) and reset dimension counter (dID)
      fID = fID + 1;
      dID = 1;
  else
      nNumFields = length(SpaceLocs);
      
      if fID == 1
         fNames{fCtr} = line( 1:SpaceLocs(1)-1 );
         fCtr = fCtr + 1;
         fDim(fCtr) = nNumFields;
      end
      
      if fID > size(D,1)
          D( end+1:end+7200, :) = zeros( 7200, size(D,2) );
      end
      
      D( fID, dID:dID+nNumFields-1 ) = sscanf( line(SpaceLocs(1)+1:end), '%f' );
      dID = dID+nNumFields;
      %nDimsPerPart(end+1) = nNumFields;
  end
end

% Make sure to close file
fclose(fid);

% Cleanup resulting data to get rid of extra rows/cols 
%   which we preallocated just in case
D = D( 1:fID, :);
keepCols = ~isnan( sum( D, 1 ) );
D = D( :, keepCols );

[basedir,~,~] = fileparts( fname );
[ColNames] = readChannelNamesFromSkeletonKey( fullfile( basedir, 'SkeletonJoints.key')  );

% ================================== POST PROCESS
% Keep only channels desired by the user, 
%  as specified by the QueryChannelNames arg

if exist( 'QueryChannelNames', 'var' )
    keepCols = [];
    for qq = 1:length( QueryChannelNames )
        needle = QueryChannelNames{qq};
        mID = find( strncmp( needle, ColNames, length(needle) ) );
        if ~isempty( mID )
            keepCols(end+1) = mID;
        else
            fprintf( 'WARNING: Did not find desired channel named %s. Skipping...\n', needle );
        end
    end
    D = D(:, keepCols );
    ColNames = ColNames( keepCols );
end

% ================================= SMOOTH ANGLE MEASUREMENTS
% Look through all channels, and correct for sharp discontinuities
%   due to angle measurements jumping past boundaries
% e.g. smooth transition from +178 deg to +182 deg could be recorded as
%                             +178 deg to -178 deg, which is awkward
didSmooth = 0;
SmoothedChannels = {};
for chID = 1:size(D,2)
    didSmoothHere = 0;
    for tt = 2:size(D, 1)
        ttDelta= D( tt, chID) - D(tt-1, chID);
        if abs(ttDelta+360)<abs(ttDelta)
            % if y_tt +360 is closer to y_tt-1 than just y_tt
            %    shift y_tt and all subsequent measurements by +360
            D(tt:end, chID) = D(tt:end, chID)+360;
            didSmoothHere= 1;
        elseif abs(ttDelta-360)<abs(ttDelta)            
            % if y_tt -360 is closer to y_tt-1 than just y_tt
            %    shift y_tt and all subsequent measurements by -360
            D(tt:end, chID) = D(tt:end, chID)-360;
            didSmoothHere= 1;            
        end
    end
    if didSmoothHere
       SmoothedChannels{end+1} = ColNames{chID};
    end
    didSmooth = didSmooth | didSmoothHere;
end
if didSmooth
    L = length( SmoothedChannels );
    MyChannels(1:2:(2*L) ) = SmoothedChannels;
    for aa = 2:2:(2*L)
        MyChannels{aa} = ', ';
    end
    SmoothSummary = strcat( MyChannels{:} );
    fprintf( 'Warning: did some smoothing on channels %s\n', SmoothSummary );
end
end % readJointAnglesAsMatrixFromAMC.m function