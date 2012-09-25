function [stateSeq, logQ_RGS] = sampleStateSeq_RGS( F, stateSeq, TS, theta, data_struct, model, objIDs, TargetPsi )

logQ_RGS = 0;
obsModel = model.obsModel;
doKeepActiveFeatsOnly = 1;
% --------------------------  loop over all unique time series objects
for ii= objIDs
    
%     if exist( 'TargetPsi', 'var' ) && ~isempty( TargetPsi )
%         % Translate from TargetPsi's feature IDs to
%         %    corresponding IDs in F, theta, prevStateSeq, etc.
%         z = TargetPsi.stateSeq(ii).z;
%         for aa = 1:length( TargetPsi.activeFeatIDs )
%             z( z == TargetPsi.activeFeatIDs(aa) ) = TargetPsi.externalFeatIDs(aa);
%         end
%         TargetPsi.stateSeq(ii).z = z;
%     end

    availFeatIDs = TS.obj(ii).availFeatIDs;    
    Kz_ii = length( availFeatIDs );
    
    pi_z = TS.obj(ii).pi_z;
    pi_z = bsxfun( @rdivide, pi_z, sum(pi_z,2)  );

    pi_init = TS.obj(ii).pi_init;
    pi_init = pi_init ./ sum( pi_init );

    T = data_struct(ii).T;    

    Kz_inds = find(F(ii,:)>0);    
    if length( Kz_inds ) == 1
        stateSeq(ii).z = double(availFeatIDs(1) ) * ones(1,T);
        continue;
    end
    
    %likelihood = compute_likelihood(data_struct(ii), theta, obsModel, Kz_inds,Kz_ii,Ks, ii);
    logLik = calcLogCondPrObsGivenTheta( data_struct(ii), theta, obsModel, Kz_inds, Kz_ii, ii, doKeepActiveFeatsOnly );
    
    normC = max( logLik, [], 1);
    logLik = bsxfun( @minus, logLik, normC );
    likelihood = exp( logLik );
    %likelihood = bsxfun(@times, likelihood, exp( normC ) );
     
    if ~exist( 'TargetPsi', 'var') || isempty( TargetPsi )
        % ------------------------------------------- Actually sample z(t)
        [js, logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, -1, randi([1 100000]) );
        z = availFeatIDs( js );
    else % ----------------------------- Calc prob of moving to Target's z
        %assert( all( find(TargetPsi.F(ii,:) ) == TargetPsi.TS.obj(ii).availFeatIDs ), 'BAD MATCH' );
        jseq = zeros(1,T);
        availTargetIDs = TargetPsi.TS.obj(ii).availFeatIDs;
        for jj = 1:length( availFeatIDs )
            jseq( TargetPsi.stateSeq(ii).z == availTargetIDs(jj) ) = jj;
        end
        [js,logQ] = SampleHMMStateSeqWithQsC( pi_z, likelihood, pi_init, jseq, randi([1 100000]) );
        z = availFeatIDs( js );
    end

    stateSeq(ii).z = z;
    
    if nargout > 1
        logQ_RGS = logQ_RGS + logQ;        
    end
   
    % Compute backwards messages:
    % pi_s = []; % Legacy parameter
    %[~, partial_marg] = MySmoothBackC( pi_z, likelihood+eps );
    %[~, partial_marg] = backwards_message_vec(likelihood, data_struct(ii).T, pi_z, pi_s);

end % loop over time series objs