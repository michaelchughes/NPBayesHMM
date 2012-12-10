% Creates the BIG motion capture dataset
%  which does NOT have ground truth action labels

DATA_DIR = '/data/liv/mhughes/data/MoCapBIG/';

MocapChannelNames = {'root.ty', 'lowerback.rx', 'lowerback.ry', 'upperneck.ry', ...
    'rhumerus.rz', 'rradius.rx','lhumerus.rz', 'lradius.rx', ...
    'rtibia.rx', 'rfoot.rx', 'ltibia.rx', 'lfoot.rx'...
    };

Preproc.obsDim = 12;
Preproc.R = 1;
Preproc.windowSize = 12;
Preproc.channelNames = MocapChannelNames;


% Create data structure
data = SeqData();

% Read in every sequence, one at a time
myFileList = dir( fullfile( DATA_DIR, 'amc', '*.amc') );
fprintf( 'Reading %d AMC files. Will take a long time...\n', length( myFileList) );
for ii = 1:length( myFileList )
    fname = myFileList(ii).name;
    
    D = readJointAnglesAsMatrixFromAMC( fullfile( DATA_DIR, 'amc', fname), Preproc.channelNames );
    fprintf( '  read in AMC file %s. # frames=%d\n', fname, size(D,1) );
    
    % Enforce zero-mean, and apply block-averaging
    D = D';
    D = bsxfun( @minus, D, mean(D,2) );
    D = preprocessMocapData(D, Preproc.windowSize );
    
    % Drop the extension ".amc" from the current sequence name string
    [~,seqName,~] = fileparts( fname );
    
    % Add the sequence to our data structure
    data = data.addSeq( D, seqName );
end
fprintf( '... read in all %d sequences\n', data.N );

saveSeqDataToPlainText( data, fullfile(DATA_DIR, 'txt') );
%fprintf( '... saved to plain text in dir %s', DATA_DIR );
