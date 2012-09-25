function [logPrF] = calcLogPrFeatureMatrix( F, gamma, c0 )
% Compute probability of observed feature assignment matrix
%   as function of the 2-parameter Indian Buffet Process 
%      using the hyperparams  gamma (mass) and c0 (conc)
% see eq. 21 of  
%   Ghahramani, Griffiths, and Sollich
%   "Bayesian nonparametric latent feature models"
%   ISBA 8th World Meeting on Bayesian Statistics, 2006


if ~exist( 'c0', 'var' )
    c0 = 1;
end

F = F(:, sum(F,1) > 0 );


    [N,K] = size( F );
        
    m = sum( F, 1 );
    
    % Find all unique column vectors in F
    [U, ~, idx_U2F] = unique( F', 'rows' );    
    nU = size(U,1 );
    logKhBang = zeros( 1, nU );
    for nn = 1:nU
        %  log ( factorial( S ) )  equiv to  gammaln( S + 1 )
        %    when S is any positive integer
       logKhBang(nn) =  gammaln(  sum( idx_U2F == nn ) + 1 );
    end
    
    % Compute log( prod_k  (N-m_k)! (m_k-1)! / N!  )    )
    %  =  sum_k log(  ( m_k - 1 )!  /  [ N (N-1) ... (N-m_k+1)  ]   )
    % factorial(n)=gamma(n+1)
    logProdFactors = zeros(1,K);
    for k = 1:K
        logProdFactors(k) = betaln( m(k), N-m(k) + c0 );
        %logProdFactors(k) = gammaln( N-m(k) + 1 ) + gammaln( m(k) ) - gammaln( N + 1 );
    end

    Harmonic_N = sum( c0./(c0 -1 + [1:N])  );
    
    logPrF = K*log( gamma ) + K*log(c0) - sum(logKhBang) - gamma*Harmonic_N  +  sum( logProdFactors );

end % calc logPr feature matrix


%     [N,K] = size( F );
%         
%     m = sum( F, 1 );
%     
%     % Find all unique column vectors in F
%     [U, ~, idx_U2F] = unique( F', 'rows' );    
%     nU = size(U,1 );
%     logKhBang = zeros( 1, nU );
%     for nn = 1:nU
%         %  log ( factorial( S ) )  equiv to  gammaln( S + 1 )
%         %    when S is any positive integer
%        logKhBang(nn) =  gammaln(  sum( idx_U2F == nn ) + 1 );
%     end
%     
%     % Compute log( prod_k  (N-m_k)! (m_k-1)! / N!  )    )
%     %  =  sum_k log(  ( m_k - 1 )!  /  [ N (N-1) ... (N-m_k+1)  ]   )
%     % factorial(n)=gamma(n+1)
%     logProdFactors = zeros(1,K);
%     for k = 1:K
%         logProdFactors(k) = gammaln( N-m(k) + 1 ) + gammaln( m(k) ) - gammaln( N + 1 );
%     end
% %     logProdFactors = zeros(1,K);
% %     for k = 1:K
% %        NUMER = factorial(  m(k) - 1   );
% %        DENOM = prod( N:-1:( N-m(k)+1 )  );
% %        logProdFactors(k) = log(NUMER) - log( DENOM );
% %     end
%     
%     Harmonic_N = sum( 1./[1:N]  );
%     
%     logPrF = K*log( gamma ) - sum(logKhBang) - gamma*Harmonic_N  +  sum( logProdFactors );
