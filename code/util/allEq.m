function status = allEq( X, Y, THR )
% Determines whether X, Y are equal to within given tolerance
%   defined on percent difference scale

if ~exist( 'THR', 'var' )
    THR = 1e-8;
end

percDiff = abs( X - Y )./abs(X);
if max(percDiff) > THR
    status = false;
else
    status = true;
end