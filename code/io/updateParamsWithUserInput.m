function pStruct = updateParamsWithUserInput(  pStruct, userInput )
% Alters default values of specified fields in pStruct struct
%   with given values provided by user
% INPUT
%   pStruct : generic structure object
%   userInput : even-length cell array of Name/Value pairs

pStruct.userInput = userInput;

% ========================================================================
%               REPLACE VALUES WITH THOSE FOUND IN userInput
% ========================================================================

assert( mod( length(userInput), 2 ) == 0, 'ModelParams Name-Value Cell Array must be EVEN!' );
for aa = 1:2:length( userInput )
    Name = userInput{aa};
    Value = force2double( userInput{aa+1} );
    
    FieldNames = textscan( Name, '%s', 'Delimiter', '.' );
    FieldNames = FieldNames{1};
    
    
    switch length( FieldNames )
        case 1
            try
                OrigVal =  pStruct.( FieldNames{1} );
                OrigDim = size( OrigVal );
                if isnumeric( OrigVal ) && ~allEq( OrigDim, size(Value) )
                    pStruct.( FieldNames{1} ) = repmat( Value, OrigDim );
                else
                    pStruct.( FieldNames{1} ) = Value;                
                end
            catch e
                pStruct.( FieldNames{1} ) = Value;
            end
            
        case 2
            try
                OrigVal =  pStruct.( FieldNames{1} ).( FieldNames{2} );
                OrigDim = size( OrigVal );
                if isnumeric( OrigVal ) && ~allEq( OrigDim, size(Value) )
                    pStruct.( FieldNames{1} ).( FieldNames{2} ) = repmat( Value, OrigDim );
                else
                    pStruct.( FieldNames{1} ).( FieldNames{2} ) = Value;
                end
            catch e
                pStruct.( FieldNames{1} ).( FieldNames{2} ) = Value;
            end
            
        case 3
            try
                OrigVal =  pStruct.( FieldNames{1} ).( FieldNames{2} ).( FieldNames{3} );
                OrigDim = size( OrigVal );
                if isnumeric( OrigVal ) && ~allEq( OrigDim, size(Value) )
                    pStruct.( FieldNames{1} ).( FieldNames{2} ).( FieldNames{3} ) = repmat( Value, OrigDim );
                else
                    pStruct.( FieldNames{1} ).( FieldNames{2} ).( FieldNames{3} ) = Value;
                end
            catch e
                pStruct.( FieldNames{1} ).( FieldNames{2} ).( FieldNames{3} ) = Value;
            end
            
    end
    
end

