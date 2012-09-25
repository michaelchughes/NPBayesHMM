function x = randdirichlet(a)
% RANDDIRICHLET   Sample from Dirichlet distribution
%    
% X = RANDDIRICHLET(A) returns a matrix, the same size as A
%       where X(:,j) is sampled from a Dirichlet(A(:,j)) distribution.
% Note: This means the *columns* of X sum to one


x = randgamma(a);
x = bsxfun( @rdivide, x, sum(x,1) );
%Z = sum(x,1);
%x = x./Z(ones(size(a,1),1),:);
