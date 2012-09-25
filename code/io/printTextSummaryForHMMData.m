function [] = printTextSummaryForHMMData( data_struct, model )

numObj = length(data_struct);
Ts = [ data_struct(:).T];

switch model.obsModel.type
    case 'Multinomial'
        data_struct(1).numVocab = length(model.obsModel.params.lambda);
        nV = data_struct(1).numVocab;
        obsDescr = sprintf( 'with %d distinct symbols\n', nV );
        if isfield( data_struct, 'nEmissions' )
            Es = [ data_struct(:).nEmissions ] ;
            obsDescr = [obsDescr sprintf('\t   between %d-%d emissions per timestep ( median = %.1f ) \n', min(Es), max(Es), median(Es)) ];
        else
            symbolInfo = sprintf('\t   all timesteps emit exactly %d symbols\n',  size( data_struct(1).obs,1 ) );
            obsDescr = [obsDescr symbolInfo];
        end
    case 'Bernoulli'
        D = size( data_struct(1).obs, 2 );
        obsDescr = sprintf( 'across %d dimensions', D );
    case 'Gaussian'
        D = size( data_struct(1).obs, 2 );
        obsDescr = sprintf( 'with %d dimensions', D );
    case 'AR-Gaussian'        
        D = size( data_struct(1).obs, 2 );
        obsDescr = sprintf( 'with %d dimensions', D );
end


fprintf( 'Dataset Summary: \n\t %d time series objects \n', numObj );
fprintf('\t   between %d-%d timesteps per object ( median = %.1f ) \n', min(Ts), max(Ts), median(Ts));
fprintf('\t   %s emissions %s\n', model.obsModel.type, obsDescr );
