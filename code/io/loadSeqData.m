function [data] = loadSeqData( datasetName, Preproc )
% Load either synthetic or real data for processing with BPHMM
%USAGE  
%   data = loadSeqData( datasetName, Preproc )
%INPUT
%   datasetName : string name
%   Preproc     : struct of preprocess params to create data
%OUTPUT
%   data  : SeqData instance

switch datasetName    
    case 'SynthAR'
        data = genToySeqData_ARGaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T, Preproc.R );
    case 'SynthGaussian'
        data = genToySeqData_Gaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T );
        
    otherwise
        error( 'Unrecognized data set name' );
end

% -------------------------------------- Preproc autoregressive data
if isfield( Preproc, 'R' ) && ~strcmp( class(data),'ARSeqData' )
    data = ARSeqData( Preproc.R, data);
end

