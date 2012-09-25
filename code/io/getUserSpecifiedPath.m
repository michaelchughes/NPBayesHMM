function saveDir = getUserSpecifiedPath( pathname )

fid = fopen( [pathname '.path'] );
saveDir = fgets( fid );
if saveDir(end)+0 == 10
    saveDir = saveDir(1:end-1); % forget about line return at the end
end
fclose( fid );