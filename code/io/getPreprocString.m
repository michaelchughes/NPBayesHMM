function [preprocStr] = getPreprocString( Preproc )
% Obtain string summary of how data was preprocessed
% Useful for creating filenames so that preproc'd data can be saved
%   and used again later without loading it from scratch   

% Traverse each field in *sorted* alphabetical order
%   This ensures similar objects always map to same string representation.
fieldNames = sort( fieldnames( Preproc ) );

preprocStr = '';
for nn = 1:length( fieldNames )
   switch( fieldNames{nn}  )
       case 'R'
           R = Preproc.( fieldNames{nn} );
           if R > 0
            preprocStr = sprintf( '%s_R%d', preprocStr, R );
           end
       case {'channelNames'}
           cIDs = Preproc.( fieldNames{nn} );
           if ~isempty( cIDs )
            preprocStr = sprintf( '%s_nCh%d', preprocStr, length(cIDs) );
           end
       case 'windowSize'
           W = Preproc.( fieldNames{nn} );
           if W > 0
            preprocStr = sprintf( '%s_W%d', preprocStr,  W);
           end
       case 'V'
           V = Preproc.( fieldNames{nn} );
           if V > 0
            preprocStr = sprintf( '%s_V%d', preprocStr,  V);
           end
       case 'timeBinSize'
           timeBinSize = Preproc.( fieldNames{nn} );
           if timeBinSize > 0
           preprocStr = sprintf( '%s_B%.2f', preprocStr, timeBinSize );
           end
       case 'nTimeBins'
           nTimeBins = Preproc.( fieldNames{nn} );
           if nTimeBins > 0
           preprocStr = sprintf( '%s_nB%d', preprocStr, nTimeBins );
           end
   end

end