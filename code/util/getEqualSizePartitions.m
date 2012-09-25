function Ns = getEqualSizePartitions( N , B )
%  Create an equitable-size partition of N objects integers into B blocks
%  Return a list of B integers, each indicating the # of objects
%    to place in each block
%EXAMPLE
%  getEqualSizePartitions( 5, 2 )
%    would return Ns=[2 3], which means a partition [1 2], [3 4 5]
%  getEqualSizePartitions( 10, 2 )
%    would return Ns=[2 2 2 2 2], which implies partition into pairs
%NOTE
%  when N not evenly divisible by B, the "extra remainder" is always added
%    back-to-front, so that Ns(end) >= Ns(end-1) >= ... Ns(1)
%    though of course, Ns(end) - Ns(1) is never more than 1

nPerB = floor(N/B);
Ns = nPerB*ones(1,B-1);
Ns(B) = N - sum( Ns );

nExcess = Ns(B) - Ns(1) - 1;

if nExcess > 0
   Ns(B-nExcess:B-1) =  Ns(B-nExcess:B-1) +1;
   Ns(B) = Ns(B) - nExcess;
end
assert( sum(Ns) == N, 'ERROR: BAD PARTITION!!' );
