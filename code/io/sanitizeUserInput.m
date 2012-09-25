function [myInput] = sanitizeUserInput( myInput )

if exist( 'myInput', 'var' ) && ~isempty( myInput )
    if ischar( myInput )
        if strmatch( '{', myInput )
            myInput = eval( myInput );
        else
            myInput = {myInput};
        end
    end
else
    myInput = {};
end