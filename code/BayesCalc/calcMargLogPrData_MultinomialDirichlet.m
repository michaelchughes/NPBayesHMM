function logPrX = calcMargLogPrData_MultinomialDirichlet( CX, lambda )

[K,V] = size( CX );
CXpL = bsxfun( @plus, CX, lambda );

%  sum_k sum_v  log Gamma( c_{k,v} + lam_v, ... c_{k,V}+lam_V  )
%  +  sum_k log Gamma( sum_v [ c_{k,v} + lam_v ]  )
logPrX = sum( sum( gammaln(  CXpL  )  ) ) ...
            - sum( gammaln( sum( CXpL,2) )  );
        
if length( lambda ) == 1
       logCLambda = K*V* gammaln( lambda ) - K*gammaln( V*lambda );
elseif numel( lambda ) == V
       logCLambda = K*sum( gammaln( lambda ) ) - K*gammaln( sum( lambda ) );
elseif size(lambda, 1) == K && size(lambda,2) == V
       logCLambda = sum(sum(gammaln(lambda))) - sum( gammaln( sum(lambda,2) ) );
end

logPrX = logPrX - logCLambda;
