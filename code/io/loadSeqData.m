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
    case {'Mocap6'}
        data = readSeqDataFromPlainText( '../data/mocap6/' );
    case {'MocapBIG'}
        data = readSeqDataFromPlainText( '/data/liv/mhughes/data/MoCapBIG/txt/' );
    otherwise
        error( 'Unrecognized data set name' );
end

% -------------------------------------- Preproc autoregressive data
% This step simply builds necessary data structs for efficient AR inference
%   including filling in XprevR data field so that for any time tt,
%       XprevR(:,tt) = [Xdata(:, tt-1); Xdata(:,tt-2); ... Xdata(:,tt-R) ]
%   allowing of course for discontinuities at sequence transitions
if isfield( Preproc, 'R' ) && ~strcmp( class(data),'ARSeqData' )
    odata = data; % keeping odata around lets debug with before/after
    data = ARSeqData( Preproc.R, odata);
end

