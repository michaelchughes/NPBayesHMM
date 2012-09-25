function [Psi, Astats, Kstats] = sampleHMMhypers( Psi, algParams )
% Sample hyperparameters (alpha,kappa) that control HMM transition distr.
% Uses a Metropolis-Hastings random walk proposal centered on current state
% INPUT:
%   Psi : input Markov Chain state
% OUTPUT
%   Psi : new Markov chain state
% DETAILS: 
%  Given hyperparams alpha (non-sticky) and kappa (sticky) params,
%   each sequence draws transition probabilities as follows:
%     eta(jj,:) ~ Gamma( alpha, alpha, ... alpha+kappa, alpha, alpha )
%     pi(jj,:) = eta(jj,:)/sum( eta(jj,:)


%------------------------ Unpack
F = Psi.F;
[N,K] = size(F);
alpha0 = Psi.TransM.prior.alpha;
kappa0 = Psi.TransM.prior.kappa;
TransM   = Psi.TransM;


% Hyperparameters for prior on kappa:
a_kappa = TransM.prior.a_kappa;
b_kappa = TransM.prior.b_kappa;
% Variance of gamma proposal:
var_kappa = algParams.HMM.var_kappa;

% Hyperparameters for prior on alpha:
a_alpha = TransM.prior.a_alpha;
b_alpha = TransM.prior.b_alpha;
% Variance of gamma proposal:
var_alpha = algParams.HMM.var_alpha;


% Build sufficient stats
%    Ki( ii ) : # total features on for object ii
%    Skk( ii ) : sum( log(   diag entries of Pz(ii)  )  )
%    Sall( ii )    : sum( log(   all entries of Pz(ii)  ) )
Ki = zeros(1,N);
Skk = zeros(1,N);
Sall = zeros(1,N);
for ii=1:N
    Ki(ii) = sum(F(ii,:));
    pi_ii = TransM.pi( ii );
    
    Skk( ii ) = sum( log( diag( pi_ii ) ) );
    Sall( ii ) = sum( log(  pi_ii(:) ) );
end

Astats.nAccept = 0;
Astats.nTotal   = algParams.HMM.Niter;

Kstats.nAccept = 0;
Kstats.nTotal   = algParams.HMM.Niter;

for nn=1:algParams.HMM.Niter
    
    %%%%%%% Sample kappa given alpha %%%%%%%
    
    % (a,b) hyperparameters of gamma prior based on fixed variance and setting
    % mean equal to previous kappa value:
    aa_kappa0 = (kappa0^2)/var_kappa;
    bb_kappa0 = kappa0/var_kappa;
    
    % Sample a proposed kappa:
    kappaP = randgamma(aa_kappa0) / bb_kappa0;
    
    % Determine log-likelihood of transition distributions given previous kappa
    % value and proposed kappa value:
    log_diff_Z = sum( ...
              Ki .* ( gammaln(alpha0*Ki+kappaP) - gammaln(alpha0*Ki+kappa0) )...
            - Ki .* ( gammaln(alpha0+kappaP) - gammaln(alpha0+kappa0)) ... 
                            + (kappaP-kappa0)*Skk ...
                    );
    % Add in prior probability of previous and proposed kappa values:
    log_diff_Prior = (a_kappa-1)*(log(kappaP)-log(kappa0))-(kappaP-kappa0)*b_kappa;
    
    % (a,b) hyperparameters of gamma prior based on fixed variance and setting
    % mean equal to proposed kappa value:
    aa_kappaP = (kappaP^2)/var_kappa;
    
    log_diff_Q  = (gammaln(aa_kappa0) - gammaln(aa_kappaP))...
                  + (aa_kappaP-aa_kappa0-1)*log(kappa0) - (aa_kappa0-aa_kappaP-1)*log(kappaP)...
                  + (aa_kappa0-aa_kappaP)*log(var_kappa);
    
    % Log accept-reject ratio:
    log_rho = log_diff_Z + log_diff_Prior + log_diff_Q;
    
    if isinf(log_rho)
        log_rho = -Inf;
    end
    rho = exp(log_rho);
    
    if rand < rho
        kappa0 = kappaP;
        Kstats.nAccept = Kstats.nAccept + 1;
    end
    % otherwise, just keep kappa0 to current value
    
    
    
    %%%%%%% Sample alpha given kappa %%%%%%%
    
    % (a,b) hyperparameters of gamma prior based on fixed variance and setting
    % mean equal to previous alpha value:
    aa_alpha0 = (alpha0^2)/var_alpha;
    bb_alpha0 = alpha0/var_alpha;
    
    % Sample a proposed alpha:
    alphaP = randgamma(aa_alpha0) / bb_alpha0;
    
    % Determine log-likelihood of transition distributions given previous alpha
    % value and proposed alpha value:
    log_diff_Z = sum( ...
                          Ki .* ( gammaln( alphaP*Ki + kappa0 ) - gammaln( alpha0*Ki + kappa0 )) ...
                        - Ki .* ( gammaln( alphaP+kappa0 )      - gammaln( alpha0+kappa0 ) ) ...
                        - Ki .* (Ki-1) .* ( gammaln( alphaP ) - gammaln( alpha0 )  ) ...
                        + ( alphaP - alpha0 ) .* Sall ...
                      );
    
    
    % Add in prior probability of previous and proposed alpha values:
    log_diff_Prior = (a_alpha-1)*(log(alphaP)-log(alpha0))-(alphaP-alpha0)*b_alpha;
    
    % (a,b) hyperparameters of gamma prior based on fixed variance and setting
    % mean equal to proposed kappa value:
    aa_alphaP = (alphaP^2)/var_alpha;
    %bb_alpha = alphaP/var_alpha; % seems to be unused
    
    log_diff_Q = (gammaln(aa_alpha0) - gammaln(aa_alphaP))...
        + (aa_alphaP-aa_alpha0-1)*log(alpha0) - (aa_alpha0-aa_alphaP-1)*log(alphaP)...
        + (aa_alpha0-aa_alphaP)*log(var_alpha);
    
    % Log accept-reject ratio:
    log_rho = log_diff_Z + log_diff_Prior + log_diff_Q;
    
    if isinf(log_rho)
        log_rho = -Inf;
    end
    rho = exp(log_rho);
    
    if rand < rho
        alpha0 = alphaP;
        Astats.nAccept = Astats.nAccept + 1;
    end
    % otherwise, just keep previous value of alpha0
    
end % end iterative loop over params

% -------------------------------  Repack
Psi.TransM.prior.alpha = alpha0;
Psi.TransM.prior.kappa = kappa0;
