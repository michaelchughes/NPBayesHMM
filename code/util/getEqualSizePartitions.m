function Ns = getEqualSizePartitions( N , B )
%  Create an equitable-size partition of N objects integers into B blocks
%  Return a list of B integers, each indicating the # of objects
%    to place in each block
%  For example, with N=5,B=2
%    we would return Ns=[2 3]

nPerB = floor(N/B);
Ns = nPerB*ones(1,B-1);
Ns(B) = N - sum( Ns );

nExcess = Ns(B) - Ns(1) - 1;

if nExcess > 0
   Ns(B-nExcess:B-1) =  Ns(B-nExcess:B-1) +1;
   Ns(B) = Ns(B) - nExcess;
end
assert( sum(Ns) == N, 'ERROR: BAD PARTITION!!' );