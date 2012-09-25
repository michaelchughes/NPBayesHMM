function [k, cs] = multinomial_single_draw( ps )
% Draw one class label from multinomial given prob vector ps
%  Basic idea is to give each class k a portion of the interval (0,1)
%   with size proportional to its probability P(Z=k)
%  We then generate a number uniformly at random within (0,1), 
%   and return the label of the bin the random number falls into

% Cumulative sum of ps give the "largest" edge of each bins
%   e.g. cumsum( [.2 .3 .5] ) = [.2  .5 1.0]
%     meaning the first class gets interval (0,.2]
%                 second      gets interval (.2, .5] etc.
% Compute cumsum(ps) > rand returns a monotonic binary vector
%   e.g.  [.2 .5 1.0] > .34567 = [0 1 1]
% The idx at which the first "1" appears from left-to-right
%   indicates the class whose bin the random value fell into
% We find this idx using the find(bs, 1, 'first') command

cs = cumsum(ps);
k = find( cs > rand*cs(end), 1);