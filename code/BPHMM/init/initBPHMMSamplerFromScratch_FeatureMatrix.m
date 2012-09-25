function F = initBPHMMSamplerFromScratch_FeatureMatrix( data_struct, initParams,outParams)
%  Features can be initialized in a few ways
%    1) unique set for each obj
%    2) unique set for each category
%    3) shared set for all objects
%           either from IBP prior
%                    or given specific number K of features to use

nObj = length( data_struct );

if isfield( initParams.F, 'nTotal' )
    K =  initParams.F.nTotal;
    F = ones( nObj, K );
    
    if outParams.doPrintHeaderInfo        
        fprintf( '\t F : %d global features shared by all objects \n', K );
    end
    
elseif isfield( initParams.F, 'nUniquePerObj' )
    Kii = initParams.F.nUniquePerObj;
    F = zeros( nObj, Kii * nObj );
    
    for ii = 1:nObj
        F(ii, (ii-1)*Kii + [1:Kii] ) = 1;
    end
    if outParams.doPrintHeaderInfo        
        fprintf( '\t F : each obj. assigned %d unique features \n', Kii );
    end
    
end