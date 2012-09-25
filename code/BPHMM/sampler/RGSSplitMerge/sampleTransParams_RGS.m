function [TS, logQ_RGS] = sampleTransParams_RGS(F, stateSeq, TS, hyperparams, model, objIDs, TargetPsi )
% Sample the transition distribution params stored in transStruct
%       state-to-state trans.   : pi_z
%       init state distribution : pi_init
% Each object ii has Pi_z matrix that is Kz_ii x Kz_ii
%   where Kz_ii := # features available for obj. ii  = sum( F(ii,:) )
% Sample HMM transition parameters under RESTRICTED settings
%   which means we
%       (1) only update params for sequences in list "objIDs"
% Executed under two circumstances:
%   (I)  actually sample from posterior
%   (II) do not sample at all, but compute prob.
%            for moving from current state to TargetPsi state

logQ_RGS = 0;

alpha0 = hyperparams.alpha;
kappa0 = hyperparams.kappa;

Zstats = getZSuffStats( F, stateSeq, model );

switch model.HMMmodel.transType
         
    case 'byObject'
        % =========================================== Sample Individual Trans Structs        
        for ii= objIDs
            
            f_ii = find( F(ii,:) );
            Kz_ii = length( f_ii );
 
            TS.obj(ii).availFeatIDs = f_ii;
            TS.obj(ii).pi_init = ones( 1, Kz_ii );

            ParamMatrix = Zstats.obj(ii).Nz     + alpha0*ones(Kz_ii,Kz_ii) + kappa0*eye( Kz_ii, Kz_ii );
            
            if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
                % Pretend we obtained current TS from target
                TS.obj(ii).pi_z = TargetPsi.TS.obj(ii).pi_z;

                Pi_z = TargetPsi.TS.obj(ii).pi_z;
                EtaSum = sum( Pi_z, 2 );
                Pi_z = bsxfun( @rdivide, Pi_z, sum(Pi_z,2)  );
            else                
                % Draw Normalized Probabilities from Dirichlet Posterior
                %   q ~ Dir( N_k + a + kappa*delta(j,k) )
                Pi_z  = randgamma(  ParamMatrix );
                Pi_z  = bsxfun( @rdivide, Pi_z, sum( Pi_z,2) );
                % Draw a scale factor for each row of Pi_z
                %   proportional to *sum* of prior parameters
                EtaSum = randgamma( (kappa0 + Kz_ii*alpha0)*ones(Kz_ii,1) );
                % Combine Dir draws with scale factor to get Gamma draws
                % via the transformation:
                %    eta_k = q_k * EtaSum  where   sum( eta_k ) = EtaSum
                TS.obj(ii).pi_z = bsxfun( @times, Pi_z, EtaSum );
            end
            
            if nargout > 1
                logQ_RGS = logQ_RGS + calcLogPrDirichlet( log(Pi_z), ParamMatrix, 1 );
                logQ_RGS = logQ_RGS + calcLogPrGamma( EtaSum, Kz_ii*alpha0 + kappa0, 1 );
            end
            
        end % loop over time series objs
        
    case 'global'
       error( 'TO DO' );
    case 'byCategory'
       error( 'TO DO (see code at end of this file for attempt 1' );
        
end % switch over trans sharing types

end % main function


% ====================================================================
function TS = initBlankTransStruct( nObj, Kz )
TS_Template = struct('pi_z',zeros(Kz,Kz, 'single'),'pi_init',zeros(1,Kz, 'single')  );
TS = repmat( TS_Template, nObj, 1 );
end