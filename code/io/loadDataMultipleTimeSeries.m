function [data] = loadDataMultipleTimeSeries( datasetName, Preproc )
% Load either synthetic or real data for processing with BPHMM
% Usage:  
%   data


switch datasetName    
    case {'SynthBern', 'SynthBernoulli', 'Bernoulli'}
        data = genSynthData_Binary( Preproc.nStates, Preproc.nDims, Preproc.nObj, Preproc.T );
    case {'SynthMultinomial', 'SynthMult', 'Mult'}
        data = genSynthData_Multinomial( Preproc.nStates, Preproc.V, Preproc.nObj, Preproc.T, Preproc.obsDim, Preproc.pEmitFavor);
    case 'SynthAR'
        data = genSynthData_ARGaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T, Preproc.R );
    case 'SynthGaussian'
        data = genSynthData_Gaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T );
        
    otherwise
        error( 'Unrecognized data set name' );
end


% -------------------------------------- Preproc autoregressive data
if isfield( Preproc, 'R' )
    data = ARSeqData( data, Preproc.R );
end

