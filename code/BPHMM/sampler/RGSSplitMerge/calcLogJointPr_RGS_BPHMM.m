function logPr = calcLogJointPr_RGS_BPHMM( Psi, data_struct, hyperparams, model, objIDs, Xstats )
% Calculate logarithm of the joint probability
%   of the state Psi output by RESTRICTED scan
% OUTPUT
%   logPr := struct with fields
%              .F
%              .z
%              .obs

logPr.eta = 0;
logPr.z = 0;
logPr.theta = 0;

featIDs = unique(Psi.activeFeatIDs);

alpha0 = hyperparams.alpha;
kappa0 = hyperparams.kappa;

logPr.F = calcLogPrFeatureMatrix( Psi.F, hyperparams.gamma, hyperparams.c );

Zstats = getZSuffStats( Psi.F, Psi.stateSeq, model, objIDs );
switch model.HMMmodel.transType % ---------------------------------------------- log p( z | F, alph, kappa)
    case 'byObject'
        for ii = objIDs
            Pz = Psi.TS.obj(ii).pi_z;
            Pi = bsxfun( @rdivide, Pz, sum(Pz,2) ); 
            logPr.eta = logPr.eta + calcLogPrGamma( Pz, [], [], alpha0, kappa0 );
            logPr.z = logPr.z + sum(sum( Zstats.obj(ii).Nz .* log(Pi) ) );
        end
    otherwise
        error( 'To Do.' );
end


if ~exist( 'Xstats', 'var' )
    Xstats = getXSuffStats( Psi.F, Psi.stateSeq, data_struct, model );
end
logPr.obs = 0;
switch model.obsModel.type
    case 'Multinomial'
        lambda = model.obsModel.params.lambda';
        logp = Psi.theta.logp;
        %logPr.theta = calcLogPrDirichlet( Psi.theta.logp( featIDs,: ), lambda, 1 );
        %logPr.obs = sum( sum( Xstats.Nkv .* logp ) );
        %Xall = vertcat(  data_struct(:).obsHist );
        %Zall = horzcat( Psi.stateSeq(:).z  );
        logPr.obs = sum( sum( Xstats.Nkv .* logp ) );
        logPr.theta = 0;        
        for kk = 1:size( Psi.F, 2 )
           if sum( Psi.F(:,kk) ) > 0
              logPr.theta = logPr.theta + calcLogPrDirichlet( logp(kk,:), lambda, 1 );
           end
        end
    case 'Gaussian'
        PP = model.obsModel.params;
        Mu = Psi.theta.Mu;
        invSigma = Psi.theta.invSigma;
        logPr.theta = 0;        
        logPr.obs = 0;
        Xall = vertcat(  data_struct(:).obs );
        Zall = horzcat( Psi.stateSeq(:).z  );
        for kk = 1:size( Psi.F, 2 )
           Xkk = Xall( Zall==kk, : );
           logPr.obs = logPr.obs + sum( calcLogPrGaussian( Xkk, Mu(kk,:), invSigma(:,:,kk) ) );
           if sum( Psi.F(:,kk) ) > 0
              logPr.theta = logPr.theta + calcLogPrNormalInvWishart( Mu(kk,:), invSigma(:,:,kk), PP );
           end
        end
    case 'AR-Gaussian'
        PP = model.obsModel.params;
        A = Psi.theta.A;
        invSigma = Psi.theta.invSigma;
        logPr.theta = 0;        
        logPr.obs = 0;
        Xall = vertcat(  data_struct(:).obs );
        XallR = vertcat(  data_struct(:).obsPrevR );
        Zall = horzcat( Psi.stateSeq(:).z  );
        for kk = 1:size( Psi.F, 2)
            keepINDS = Zall==kk;
            Xkk = Xall(   keepINDS, : );
            XkkR = XallR( keepINDS, : );
            logPr.obs = logPr.obs + sum( calcLogPrGaussian( Xkk, XkkR*A(:,:,kk)', invSigma(:,:,kk)  ) );

            if sum( Psi.F(:,kk) ) > 0
                logPr.theta = logPr.theta + calcLogPrMatrixNormalInvWishart( A(:,:, kk), invSigma(:,:,kk), PP );
            end
        end
    otherwise
        error( 'to do.' );
end


logPr.all = logPr.obs + logPr.z + logPr.F + logPr.eta + logPr.theta;

end