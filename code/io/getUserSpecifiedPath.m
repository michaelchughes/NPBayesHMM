function saveDir = getUserSpecifiedPath( pathname )

if ~exist( [pathname '.path'], 'file' )
   error( sprintf('ERROR: required path definition %s.path does not exist. Please create it.', pathname ) );
end

fid = fopen( [pathname '.path'] );
saveDir = fgets( fid );
if saveDir(end)+0 == 10
    saveDir = saveDir(1:end-1); % forget about line return at the end
end
fclose( fid );

if ~exist( saveDir, 'dir' )
    [~,~] = mkdir( saveDir );
end