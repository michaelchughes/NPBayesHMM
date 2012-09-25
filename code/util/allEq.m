function status = allEq( X, Y, THR )
% allEq() : rapidly assert X and Y are (almost exactly) the same
% Check if matrices X, Y are equal to within given tolerance
%   defined on the percent difference scale
%INPUT: 
%   X, Y: N x D matrices
%   THR: numeric tolerance for equality (default 1e-8)
%OUTPUT
%  status : boolean
%     true if EVERY entry has perc diff <= THR
%     false o.w.

if ~exist( 'THR', 'var' )
    THR = 1e-8;
end

percDiff = abs( X - Y )./abs(X);
if max(percDiff) > THR
    status = false;
else
    status = true;
end
