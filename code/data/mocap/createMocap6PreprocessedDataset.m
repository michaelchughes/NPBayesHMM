% Creates the 6 sequence motion capture dataset
%  from subjects 13 and 14 of the CMU database


DATA_DIR = '../data/mocap6/';

MocapChannelNames = {'root.ty', 'lowerback.rx', 'lowerback.ry', 'upperneck.ry', ...
    'rhumerus.rz', 'rradius.rx','lhumerus.rz', 'lradius.rx', ...
    'rtibia.rx', 'rfoot.rx', 'ltibia.rx', 'lfoot.rx'...
    };

Preproc.nObj = 6;
Preproc.obsDim = 12;
Preproc.R = 1;
Preproc.windowSize = 12;
Preproc.channelNames = MocapChannelNames;

fprintf( 'Reading %d AMC files. Will take a long time...\n', length( myFileList) );
GT = load( fullfile( DATA_DIR, 'mat', 'ExerciseGroundTruth.mat' ) );
GT = GT.GT;
myFileList = dir( fullfile( DATA_DIR, 'raw', 'amc', '*.amc') );

% Create data structure
data = SeqData();

% Read in every sequence, one at a time
fprintf( 'Reading %d AMC files. Will take a long time...\n', length( myFileList) );
for ii = 1:length( myFileList )
    fname = myFileList(ii).name;
    assert( strcmp( fname, GT(ii).fname ), 'ERROR: Not matched to ground truth!' );
    
    D = readJointAnglesAsMatrixFromAMC( fullfile( DATA_DIR, 'raw', 'amc', fname), Preproc.channelNames );
    fprintf( '  read in AMC file %s. # frames=%d\n', fname, size(D,1) );
    
    % Enforce zero-mean, and apply block-averaging
    D = D';
    D = bsxfun( @minus, D, mean(D,2) );
    D = preprocessMocapData(D, Preproc.windowSize );
    
    % Drop the extension ".amc" from the current sequence name string
    [~,seqName,~] = fileparts( fname );
    
    % Add the sequence to our data structure
    data = data.addSeq( D, seqName,  GT(ii).true_labels );
end
fprintf( '... read in all %d sequences\n', data.N );

saveSeqDataToPlainText( data, DATA_DIR );
fprintf( '... saved to plain text in dir %s', DATA_DIR );
