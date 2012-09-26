function [CNames] = readChannelNamesFromSkeletonKey( fname )

% Open file
fid=fopen(fname, 'rt');
if (fid == -1)
  error('ERROR: Cannot open file %s.\n', fname);
end;

% Read lines until we've skipped past the header 
%   assuming it ends with the line ":DEGREES"
line=fgetl(fid);
if strcmp(line(1),'#')
  while strcmp(line(1),'#')
      line=fgetl(fid);
  end
end

CNames = {};
doReadFreshLine = 0;
while ~feof(fid)

  if doReadFreshLine
      line = fgetl( fid );
  end
  doReadFreshLine = 1;
  
  fields = regexp(line,' ','split');  
  
  PartName = fields{1};
  
  for aa = 2:length( fields )
      CNames{end+1} = [PartName '.' fields{aa}];
  end
  

end

