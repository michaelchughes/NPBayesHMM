function [data] = readSeqDataFromPlainText( outputPath )
% Read a collection of sequence data from plain text
%INPUT:
%  data : SeqData object (or any subclass)
%  outputPath : string valid path to save results to


if ~exist( outputPath, 'dir' )
   error( 'ERROR: no such directory' ); 
end


namepath = fullfile(outputPath, 'SeqNames.txt' );
namefid = fopen(namepath);
Names = textscan( namefid, '%s' );
Names = Names{1};
fclose(namefid);

ztruepath = fullfile(outputPath, 'zTrue.dat' );
if exist( ztruepath, 'file' )
    fid = fopen( ztruepath );
    line = fgetl(fid);
    ii=1;
    while line >= 0 %true until EOF reached
       zTrue{ii} = str2num( line );
       line = fgetl(fid); 
       ii=ii+1;       
    end
    fclose( fid);
end
dirList = dir( fullfile(outputPath, '*.dat') );

fileNames = {dirList(:).name};

data = SeqData();
for ii = 1:length( Names )
    queryName = [Names{ii} '.dat'];
    
    aa = strmatch( queryName, fileNames );
    
    if isempty( aa )
        error([ 'ERROR: Cannot find data for sequence named ' queryName] );
    end
    
    Xseq = importdata( fullfile(outputPath, queryName) );
    
    if exist( 'zTrue', 'var' )
        data=data.addSeq( Xseq', Names{ii}, zTrue{ii} );
    else
        data=data.addSeq( Xseq', Names{ii} );
    end
end