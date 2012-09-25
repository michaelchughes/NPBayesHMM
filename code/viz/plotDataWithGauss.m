function [] = plotDataWithGauss( X, z, MyColors )

if ~exist( 'MyColors', 'var' );
    MyColors = repmat( [0 0 0], max( unique(z) ), 1 );
end;

z = z(:)';
for kk = unique(z)
   INDS = (z == kk);
   plot( X(INDS,1), X(INDS,2), '+', 'MarkerSize', 4, 'Color', MyColors(kk,:) );
end