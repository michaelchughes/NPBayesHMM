function [Psi, Stats] = sampleIBPhypers(Psi, algParams)
% Resample two hyperparams of the Indian Buffet Process (mass & conc)
%   given the observed feature matrix F and Gamma priors on these params
% Procedes by conditional updates to one given the other
%   e.g. gamma ~ p( gamma | F, c0    )
%         c0   ~ p( c0    | F, gamma )
% Sampling Details
%   gamma0 (mass param )
%       under a Gamma( a_g, b_g ) prior
%       has a closed-form Gamma posterior distribution
%         with parameters   a = a_g  + K 
%                           b = b_g  + sum_i=1^N  c0/(c0 + i - 1 )
%   c0 (conc param)
%       under a Gamma( a_c, b_c ) prior
%       can be updated via a Met-Hastings update 
%         with proposal  q( c* | c ) = Gamma( mean=c, var = var_c )
%           where var_c is a param. of the sampling algorithm
%            and is defined within the algParams struct
% References:
%   Gharamani, Griffiths, and Solich. Esp. eq. 21


% =============================== UNPACK
F = Psi.F;
gamma0 = Psi.bpM.gamma;
a_gamma = Psi.bpM.prior.a_mass;
b_gamma = Psi.bpM.prior.b_mass;

C0     = Psi.bpM.c;
a_c0 = Psi.bpM.prior.a_conc;
b_c0 = Psi.bpM.prior.b_conc;

var_c0 = algParams.BP.var_c;

nObj = size(F,1);
Kplus = sum(sum(F,1)>0);
featCounts = sum(F,1);

Stats.nAccept = 0;
Stats.nTotal   = algParams.BP.Niter;

for iter = 1:algParams.BP.Niter
    
    if algParams.BP.doSampleMass
        % ------------------------------------------------ Update gamma0 | c0
        Harmonic_Cur = sum( C0 ./ ( C0 - 1 + (1:nObj) )  );

        a_POST = a_gamma + Kplus;
        b_POST = b_gamma + Harmonic_Cur;
        
        % Draw X ~ Gamma(a,b) is equiv. in distribution to
        %      Y/b  when Y ~ Gamma(a,1)
        gamma0 = randgamma( a_POST ) / b_POST;
    end
    
    % ------------------------------------------------ Update c0 | gamma0
    if algParams.BP.doSampleConc
        
        a_q = (C0^2) / var_c0;
        b_q = C0   / var_c0;
        
        propC0 = randgamma( a_q ) / b_q;
        Harmonic_New = sum(  propC0 ./ ( propC0 - 1 + [1:nObj] ) );
        
        logProd_New = sum( betaln( featCounts, nObj - featCounts + propC0  ) );
        logProd_Cur = sum( betaln( featCounts, nObj - featCounts + C0  )  );
        
        logPrF_New = Kplus*log( propC0 ) - gamma0*Harmonic_New + logProd_New;
        logPrF_Cur = Kplus*log( C0 )     - gamma0*Harmonic_Cur + logProd_Cur;
        
        logPrC0_New = (a_c0-1)*log( propC0 ) - b_c0*propC0;
        logPrC0_Cur = (a_c0-1)*log( C0 )     - b_c0*C0;
        
        a_qR = propC0^2/var_c0;
        b_qR = propC0/var_c0;
        
        logQ = a_q*log( b_q ) - gammaln( a_q ) + (a_q-1)*log( propC0 ) - b_q*propC0;        
        logQRev = a_qR*log( b_qR ) - gammaln( a_qR ) + (a_qR-1)*log( C0 ) - b_qR*C0;
        
        logPrAccept = logPrF_New - logPrF_Cur + logPrC0_New - logPrC0_Cur + logQRev - logQ;
        
        rho = exp( logPrAccept );
        if rand < rho
            C0 = propC0;            
            Stats.nAccept = Stats.nAccept + 1;
        end
        
    end

end


% ========================== REPACK 
Psi.bpM.gamma = gamma0;
Psi.bpM.c = C0;

end % main function

