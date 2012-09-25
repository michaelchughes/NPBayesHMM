function d = force2double( s , doLeaveAlone)
% Coerce input s to be of type double
%   Can handle strings or any numeric input 

if ~exist( 'doLeaveAlone', 'var' )
    doLeaveAlone = 0;
end

if ischar( s )
    if doLeaveAlone
        d = s;
        return;
    else
        
        if s(1) == '[' && s(end) == ']'
            d = eval( s );
        elseif s(1) == '{' && s(end) == '}'
            s = [ '[' s(2:end-1) ']' ];
            d = eval( s );
        elseif strfind( s, '/' );
            % process fraction inputs, eg. '1/25' as 0.04, etc...
            d = eval( s );  
        else
            d = str2double(s);
        end
    end
    
    if isnan( d )
        d = s;
    end
    
elseif strcmp( class(s), 'double' )
    d = s;
elseif isnumeric(s)
    d = str2double( num2str(s) );
elseif strcmp( class(s), 'function_handle'  )
    d = s;
else
    if doLeaveAlone
        d = s;
    else
        error( 'Unsupported type!');
    end
end