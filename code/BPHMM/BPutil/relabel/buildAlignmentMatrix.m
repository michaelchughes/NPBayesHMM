function [A uTrue uEst] = buildAlignmentMatrix( zTrue, zEst )
% Given true and estimated label sequence, construct align matrix A
%   where A( ke, kt ) = fraction of labels ke aligned 
%                          from ESTIMATED label ke to TRUE label kt
% Note that each ROW of A will sum to 1

uTrue = unique( zTrue );
uEst  = unique( zEst  );

for ue = 1:length( uEst )
    ts = zEst == uEst(ue);
    A(ue,:) = histc( zTrue(ts),  uTrue );
end
A = bsxfun( @rdivide, A, sum(A,2)  );