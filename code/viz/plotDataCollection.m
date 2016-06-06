function [] = plotDataCollection( data )

N = min( 10, data.N );
ks = unique(data.zTrueAll );
%Colors = jet(length(ks));
Colors = get(0,'defaultAxesColorOrder');
Colors( 8 , : ) = [ 1 0.5 0];
for ii = 1:N
   zTrue = data.zTrue(ii);
   
   X = data.seq(ii);
   
   ts = 1:data.Ts(ii);
   
   subplot( N, 1, ii );
   hold all;
   for kk = 1:length(ks)
      kkINDS = zTrue == ks(kk);
      if ~isempty(kkINDS)
      plot( ts(kkINDS), X(1,kkINDS), '.', 'Color', Colors(kk,:), 'MarkerSize', 20 );
      plot( ts(kkINDS), X(2,kkINDS), '.', 'Color', Colors(kk,:), 'MarkerSize', 12 );
      end
   end
   ylim( [-2 2] );
   set( gca, 'YTick', [] );
   set( gca, 'FontSize', 20 );
end


end