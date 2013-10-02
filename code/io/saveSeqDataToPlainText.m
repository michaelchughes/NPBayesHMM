function [] = saveSeqDataToPlainText( data, outputPath )
% Save a collection of sequence data to plain text
%INPUT:
%  data : SeqData object (or any subclass)
%  outputPath : string valid path to save results to


[~,~] = mkdir( outputPath );
for ii = 1:data.N
   
    seqName =  data.name(ii);
    
    fid = fopen( fullfile( outputPath, [seqName '.dat'] ) , 'w');
    
    FMT_STR = [repmat('%f ', 1, data.D) '\n'];
    
    X = data.seq(ii);
    for tt = 1:data.Ts(ii)
        fprintf( fid, FMT_STR, X(:,tt) );
    end
    
    fclose(fid);    
end

if isfield( data, 'zTrue' ) && ~isempty( data.zTrue(1) )
   fid = fopen( fullfile( outputPath, ['zTrue.dat'] ), 'w');
   namefid = fopen( fullfile( outputPath, ['SeqNames.txt'] ), 'w');
   for ii = 1:data.N
       zTrue = data.zTrue(ii);
       FMT_STR = [repmat('%d ', 1, data.Ts(ii)) '\n'];
       
       fprintf( fid, FMT_STR, zTrue);
       
       fprintf( namefid, '%s\n', data.name(ii) );       
   end
   fclose(fid);
   fclose( namefid );
else
    % Save plain text file of sequence names (to define specific ordering)
    namefid = fopen( fullfile( outputPath, ['SeqNames.txt'] ), 'w');
    for ii = 1:data.N       
        fprintf( namefid, '%s\n', data.name(ii) );
    end
    fclose( namefid );
end