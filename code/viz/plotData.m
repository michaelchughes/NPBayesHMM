function [] = plotData( data, ii )

T = data.Ts(ii);
M = 10;

xs = 1:T;
ys = linspace( -3, 3, M);

hold all;
hIM = imagesc( xs, ys, repmat(data.zTrue(ii), M, 1), [1 max( data.zTrueAll)] );
set( hIM, 'AlphaData', 0.65 );
X = data.seq(ii);
plot( xs, X(1,:), 'k.-' );
plot( xs, X(2,:), 'r.-' );
title( ['Sequence ' num2str(ii)], 'FontSize', 20 );

axis( [1 T ys(1) ys(end)] );

end