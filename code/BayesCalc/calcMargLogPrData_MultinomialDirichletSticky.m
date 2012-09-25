function logPrX = calcMargLogPrData_MultinomialDirichletSticky( CX, alpha, kappa )

[K,V] = size(CX );
lambda = alpha*ones(K,K) + kappa*eye(K);
logPrX = calcMargLogPrData_MultinomialDirichlet( CX, lambda );

end