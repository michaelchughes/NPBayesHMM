function [TS] = getEta_PriorMean(F, stateCounts, hyperparams, model, allObjIDs, featIDs )
% Sample HMM transition parameters under RESTRICTED settings
%   which means we
%       (1) only update params for sequences in list "objIDs"
% Executed under two circumstances:
%   (I)  actually sample from posterior
%   (II) do not sample at all, but compute prob.
%            for moving from current state to TargetPsi state
% OUTPUT
%    Under mode (I),
%            TS.obj(ii).availFeatIDs <-- find( F(ii,:) )
%            TS.obj(ii).pi_z  <-- Posterior(pi_z | stateCounts, priors )
%    Under mode (II),
%            TS.obj(ii)  <-- set to TargetPsi.TS.obj(ii) always

alpha0 = hyperparams.alpha;
kappa0 = hyperparams.kappa;

switch model.HMMmodel.transType
    case 'byObject'
        % =========================================== Sample Individual Trans Struct
        TSperObj = struct( 'availFeatIDs', [], 'pi_init', [], 'pi_z', [] );
        TS.obj = repmat( TSperObj, size(F,1), 1 );
        for ii=allObjIDs
            % Obtain current avail feat IDs for this object
            f_ii = sort( union(find( F(ii,:) ), featIDs )  );

            TS.obj(ii).availFeatIDs = f_ii;
            Kz_ii = length( f_ii );

            TS.obj(ii).pi_init = ones( 1, Kz_ii );
            
            TS.obj(ii).pi_z = alpha0 +  kappa0*eye( Kz_ii, Kz_ii );
            
        end % loop over time series objs
 
    case 'byCategory'
        % ===========================================  Sample By Class Trans Struct
        error( 'TO DO' );
    otherwise
        error( 'TO DO' );
        
end % switch over trans sharing types


end % main function