function [data_struct] = loadDataForBPHMM( datasetName, splitCat, Preproc, settings )
% Load either synthetic or real data for processing with BPHMM
% Usage:  
%   data_struct = loadDataForBPHMM( <datasetName>, <splitName>, ...
%                                     dataInfo, settings )
%  Should be in struct array format with following fields
%   data_struct(ii)
%       .T
%       .nEmissions
%       .obs
%       .true_labels


doLoadFromFile = 1;

switch datasetName
    case 'MoCapBIG'
        
        DATA_DIR = '/data/liv/mhughes/data/MoCapBIG/';        
        addpath( genpath( '~/MoCap/' ) );
        
        matfilename = sprintf( 'MoCapSensorData%s.mat', getPreprocString( Preproc ) );
        matfilepath = fullfile( DATA_DIR, 'MAT', matfilename );
        if exist( matfilepath, 'file' )
            MAT = load(  matfilepath );
            data_struct = MAT.data_struct;
            fprintf( '... read AMC data from MAT file: %s\n', matfilepath );

        else
            myFileList = dir( fullfile( DATA_DIR, 'amc', '*.amc') );
            for ii = 1:length( myFileList )
                % replace 'MoCap/data/' path with '' (nothing)
                fname = myFileList(ii).name;
                D = readAMC( fullfile( DATA_DIR, 'amc', fname), 1 );
                fprintf( '  read in AMC file %s. # frames=%d\n', fname, size(D,1) );
                % Transpose to apply same preprocessing as in orig. BP-HMM
                D = D';
                if isfield( Preproc, 'channelIDs' )
                    D = D( Preproc.channelIDs, : );
                end;
                meanD = mean(D,2);
                D = bsxfun( @minus, D, meanD );
                D = preprocessData(D, Preproc.windowSize );
                
                % Transpose back, since we want it in T x obsDim format
                data_struct(ii).obs = D';
                data_struct(ii).fileName = fname;
                data_struct(ii).T = size( data_struct(ii).obs, 1);
            end
            save(  matfilepath , 'data_struct' );
            fprintf( '... wrote to file for future quick loading: %s\n', matfilepath );
        end
        doLoadFromFile = 0;
    case 'MoCap6'
        if ismac
            DATA_DIR = '~/freeware/matlab/BPARHMMtoolbox/data/';
            Q = load( '~/freeware/matlab/BPARHMMtoolbox/Erik/exerciseRoutines.mat' );
            addpath( genpath( '~/freeware/matlab/BPARHMMtoolbox/Erik/' ) );
        else
            DATA_DIR = '/data/liv/mhughes/data/MoCap6/';        
            Q = load( fullfile( DATA_DIR, 'MAT', 'exerciseRoutines.mat' ) );
            addpath( genpath( '~/MoCap/' ) );
        end
                
        matfilename = sprintf( 'MoCapSensorData%s.mat', getPreprocString( Preproc ) );
        matfilepath = fullfile( DATA_DIR, 'MAT', matfilename );
        if exist( matfilepath, 'file' )
            MAT = load(  matfilepath );
            data_struct = MAT.data_struct;
            fprintf( '... read AMC data from MAT file: %s\n', matfilepath );
        else
            
            RawData = Q.Data;
            for ii = 1:length( RawData )
                % replace 'MoCap/data/' path with '' (nothing)
                fname = strrep( Q.fname{ii}, '../MoCap/data/', '' );
                if ismac
                    D = readAMC( fullfile( DATA_DIR, fname), 1 );
                else
                    D = readAMC( fullfile( DATA_DIR, 'amc', fname), 1 );
                end
                fprintf( '  read in AMC file %s. # frames=%d\n', fname, size(D,1) );
                % Transpose to apply same preprocessing as in orig. BP-HMM
                D = D';
                if isfield( Preproc, 'channelIDs' )
                    D = D( Preproc.channelIDs, : );
                end;
                meanD = mean(D,2);
                D = bsxfun( @minus, D, meanD );
                D = preprocessData(D, Preproc.windowSize );
                
                % Transpose back, since we want it in T x obsDim format
                data_struct(ii).obs = D';
                data_struct(ii).fileName = fname;
                data_struct(ii).T = size( data_struct(ii).obs, 1);
                data_struct(ii).true_labels = Q.true_labels{ii};
                
            end 
            data_struct(3).true_labels( data_struct(3).true_labels==10 ) = 12;

            save(  matfilepath , 'data_struct' );
            fprintf( '... wrote to file for future quick loading: %s\n', matfilepath );
        end
        doLoadFromFile = 0;

    case 'ChromatinStates'
        X = load( '/data/liv/mhughes/data/ChromatinStates/HMMdata/MarkerSeqData_custom.mat' );
        data_struct = X.data_struct;
        doLoadFromFile = 0;
    case 'HMMSynthEasy'
        data_struct = genHMMSynthData_Multinomial( Preproc.nStates, Preproc.V, Preproc.nObj, Preproc.T, Preproc.obsDim, Preproc.pEmitFavor);
        doLoadFromFile = 0;
    case 'SynthBernoulli'
        data_struct = genSynthData_Binary( Preproc.nStates, Preproc.nDims, Preproc.nObj, Preproc.T );
        doLoadFromFile = 0;
    case 'Synth'
        if Preproc.nObj == 0
            Preproc.nObj = 50;
        end
        data_struct = genSynthData_Multinomial( Preproc.nStates, Preproc.V, Preproc.nObj, Preproc.T, Preproc.obsDim, Preproc.pEmitFavor);
        doLoadFromFile = 0;
    case 'SynthAR'
        data_struct = genSynthData_ARGaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T, Preproc.R );
        doLoadFromFile = 0;
    case 'SynthGaussian'
        if Preproc.nObj == 0
            Preproc.nObj = 25;
        end
        data_struct = genSynthData_Gaussian( Preproc.nStates, Preproc.obsDim, Preproc.nObj, Preproc.T );
        doLoadFromFile = 0;
    case 'SynthNoisy'
        data_struct = genSynthData_Multinomial_WithObjNoise( Preproc.nStates, Preproc.V, Preproc.nObj, Preproc.T, Preproc.obsDim, Preproc.pEmitFavor, Preproc.rNoise);
        doLoadFromFile = 0;
        
    case 'SynthHBP'
        data_struct = genSynthData_HBP( Preproc.nStates, Preproc.V, Preproc.nObjPerCat, Preproc.nCat, Preproc.pEmitFavor, Preproc.obsDim, Preproc.T);
        doLoadFromFile = 0;
    case 'SynthHBP4'
        DATA_DIR = '/data/liv/mhughes/data/SynthHBP4/';


    case 'SuperSynth'
        
        DATA_DIR = '/data/liv/mhughes/data/SuperSynth/';
        
    case 'SuperSynthC2'
        DATA_DIR = '/data/liv/mhughes/data/SuperSynthC2/';
    case 'KTHss2'
        DATA_DIR = '/data/liv/mhughes/data/KTHss2/';
    case 'KTHsubseq'
        DATA_DIR = '/data/liv/mhughes/data/KTHsubseq/';
        
    case 'KTHremix'
        DATA_DIR = '/data/liv/mhughes/data/KTHremix/';
    case 'OlympicSports'
        DATA_DIR = '/data/liv/mhughes/data/OlympicSports/';
    otherwise
        try
            DATA_DIR = ['/data/liv/mhughes/data/' datasetName];
            GT = loadGroundTruth( DATA_DIR );
        catch exception
            error( 'Unrecognized data set name' );
        end
end

if doLoadFromFile
    if exist( 'Preproc','var' ) && ~isempty( fieldnames(Preproc) )
        filename = sprintf( 'TimeSeriesBoF%s.mat', getPreprocString( Preproc )  );
        datafilepath = fullfile( DATA_DIR, 'TimeSeriesBoF', Preproc.detType, Preproc.featureName, splitCat, filename );
    else
        datafilepath = fullfile( DATA_DIR, 'TimeSeriesBoF', splitCat, 'TimeSeriesBoF.mat');
    end
    
    if exist( 'settings', 'var' ) && isfield( settings, 'Data' ) && isfield( settings.Data, 'jobID' ) && isfield( settings.Data, 'taskID' )
        CorpusData.BoFseq = loadSamplerInfo( settings.Data.jobID, settings.Data.taskID, '', {'data_struct'} );
    elseif ~exist( datafilepath, 'file' );
        data_struct = -1;
        return;
    else
        CorpusData = load( datafilepath );    
    end
    
    GT = loadGroundTruth( DATA_DIR );
    
    % Add category label info
    if ~isfield( CorpusData.BoFseq, 'classID' )
        if isfield( GT, 'doLabelVideo') 
            if GT.doLabelVideo
                for ii = 1:length( CorpusData.BoFseq )
                    classID = GT.y.( splitCat )( ii , : );
                    CorpusData.BoFseq(ii).classID = classID;
                end
            end
        else
            for ii = 1:length( CorpusData.BoFseq )
                classID = GT.y.( splitCat )( ii , : );
                CorpusData.BoFseq(ii).classID = classID;
            end
        end
    end
    
    data_struct = struct();
    % If necessary, filter out un-requested Categories
    if exist( 'Preproc', 'var')
        if isfield( Preproc, 'doConcatSeq' ) && Preproc.doConcatSeq > 0
            concatData = struct();
            dd = 0;
            for cc = 1:GT.nCategories
                if isfield( Preproc, 'Categories' )
                   if isempty( strmatch( GT.CategoryNames{cc}, Preproc.Categories )  )
                       continue;
                   end
                end
                data_struct = CorpusData.BoFseq;
                % Find all objects with this category label
                ccObjs = find( GT.y.(splitCat) == cc )';
                concatObs = {};
                nEmissions = [];
                for ii = ccObjs
                    T_ii = length( data_struct(ii).obs );
                    for tt = 1:T_ii
                        concatObs{ end+1 } = data_struct(ii).obs{tt};
                        nEmissions( end+1 ) = length( data_struct(ii).obs{tt} );
                    end
                end
                dd = dd + 1;
                concatData(dd).nEmissions = nEmissions;
                concatData(dd).obs = concatObs;
                concatData(dd).T = length( concatObs );
            end
            data_struct = concatData;
            
        elseif isfield( Preproc, 'Categories' ) && iscell( Preproc.Categories )
            for cc = 1:length( Preproc.Categories )
                catName = Preproc.Categories{cc};
                
                catID = strmatch( catName, GT.CategoryNames );
                catIDs(cc) = catID;
                if GT.doExclusive
                    %objIDs.( catName ) = find(  GT.y.(splitCat) == catID  );
                    objIDs.( catName ) = find( [CorpusData.BoFseq(:).classID] == catID  );
                else
                    % TO DO: fill in this for non-exclusive data
                end
            end
            
            % Sort by catID, 
            %    so we always process same categories in same order
            % e.g. if user A asks for {'walking', 'running'}
            %       and user B asks for {'running', 'walking'}
            %  they will get EXACT same dataset in ordering of data items
            [~, sortIDs] = sort( catIDs );
            
            for cc = sortIDs
                catName = Preproc.Categories{ cc  };

                if Preproc.nObj < 0
                    permStream = RandStream.create( 'mt19937ar', 'Seed', sum( catName )  );
                    permIDs = randperm(permStream,  length( objIDs.(catName)  )  );
                    LL = -Preproc.nObj;                    
                else
                    if Preproc.nObj == 0
                        LL = length( objIDs.(catName) );
                    else
                        LL = Preproc.nObj;
                    end
                    permIDs = 1:length( objIDs.(catName)  );
                end
                LL = min( LL, length( objIDs.(catName)  )  );                    
                objIDs.( catName ) = objIDs.( catName )( permIDs(1:LL) );                    
                
                
                for aa = 1:length( objIDs.( catName )  )
                    bb = objIDs.( catName )( aa );
                    if isempty( fieldnames( data_struct )  )
                        data_struct = CorpusData.BoFseq(bb);
                    else
                        data_struct( end + 1) = CorpusData.BoFseq( bb );
                    end                    
                end
                
            end
            Lprev = 0;
            for cc = sortIDs
               catName = Preproc.Categories{ cc  };
               for aa = 1:length( objIDs.( catName )  )
                    bb = objIDs.( catName )( aa );
                    data_struct( Lprev + aa ).objID = bb;                
               end
               Lprev = Lprev + length( objIDs.( catName ) );
            end
        elseif isfield( Preproc, 'nObj' ) && Preproc.nObj > 0
            LL = min( Preproc.nObj, length( CorpusData.BoFseq  )  );
            permStream = RandStream.create( 'mt19937ar', 'Seed', 918 );
            permIDs = randperm(permStream,  length( CorpusData.BoFseq  )  );

            data_struct = CorpusData.BoFseq( permIDs(1:LL)  );
            for ii = 1:length( data_struct )
                data_struct(ii).objID = permIDs(ii );
            end
        else
            data_struct = CorpusData.BoFseq;
        end
    end
    
    data_struct(1).DATA_DIR = DATA_DIR;    
    Preproc.DATA_DIR = DATA_DIR;
    
end

% -------------------------------------- Compute observation histogram at each time point
if isfield( Preproc, 'V' )
    for ii = 1:length( data_struct )
        obsHist = zeros( data_struct(ii).T,   Preproc.V );
        for tt = 1:data_struct(ii).T
            if iscell( data_struct(ii).obs )
                obsHist(tt, :) = histc( data_struct(ii).obs{tt},  1:Preproc.V );
            else
                obsHist(tt, :) = histc( data_struct(ii).obs(:,tt),  1:Preproc.V );
            end
        end
        data_struct(ii).obsHist = sparse( obsHist );
    end
elseif isfield( Preproc, 'R' )
    R = Preproc.R;
    for ii = 1:length( data_struct )
        [T D] = size( data_struct(ii).obs );
        obsPrevR = zeros( T-1,  D* R );
        xprev = zeros( 1, D*R );
        for rr = 1:R
            xprev( (R-rr)*D + (1:D) ) = data_struct(ii).obs(rr,:);
        end
        % Omit first R observations from actual dataseries,
        %   but store them in the "Previous" data xprev
        for tt = R+1:data_struct(ii).T
            obsPrevR(tt-R,:) = xprev;
            xprev = [data_struct(ii).obs(tt,:) xprev( 1:(end-D) ) ];
        end
        data_struct(ii).T = data_struct(ii).T - R;
        data_struct(ii).obsPrevR = obsPrevR( 1:data_struct(ii).T, : );
        data_struct(ii).obs      = data_struct(ii).obs( R+1:end, : );
        if isfield( data_struct(ii), 'true_labels' )
            data_struct(ii).true_labels = data_struct(ii).true_labels(R+1:end );
        end
    end
end


% ----------------------------------------------------- Init fixed z seq.
% if exist( 'settings', 'var') && isfield( settings, 'Init') && settings.Init.nUniqueF > 0
%     nBlocks = max( settings.Init.nUniqueF, 1);
%     
%     z_max = 0;
%     for ii = 1:length(data_struct)
%         
%         % Form initial mode sequences to simply block partition each
%         % time series into 'Ninit' features.  Time series are given
%         % non-overlapping feature labels:
%         T = data_struct(ii).T;
%         
%         init_blocksize = floor(T/nBlocks);
%         z_init = [];
%         for i=1:nBlocks
%             z_init = [z_init i*ones(1,init_blocksize)];
%         end
%         z_init(nBlocks*init_blocksize+1:T) = nBlocks;
%         data_struct(ii).z_init = z_init + z_max;
%         
%         z_max = max(data_struct(ii).z_init);
%     end
% end
