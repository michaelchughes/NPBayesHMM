function logPr = calcLogJointPr_RestrictedSet_BPHMM( Psi, data_struct, hyperparams, model, objIDs, Xstats )
% Calculate logarithm of the joint probability
%   of the state Psi output by RESTRICTED scan
% OUTPUT
%   logPr := struct with fields
%              .F
%              .z
%              .obs


logPr.F = calcLogPrFeatureMatrix( Psi.F, hyperparams.gamma, hyperparams.c );

logPr.z = 0;
%[stateCounts, INDS] = getStateSeqSuffStats( Psi.stateSeq, Psi.F, model );
Zstats = getZSuffStats( Psi.F, Psi.stateSeq, model, objIDs );
switch model.HMMmodel.transType % ---------------------------------------------- log p( z | F, alph, kappa)
    case 'byObject'
        for ii = objIDs
            logPr.z = logPr.z + calcMargLogPrData_MultinomialDirichlet(  Zstats.obj(ii).Nz, [], 1, hyperparams.alpha, hyperparams.kappa );
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
        logPr.obs = calcMargLogPrData_MultinomialDirichlet( Xstats.Nkv, model.obsModel.params.lambda', 1 );
    case 'Bernoulli'
        error( 'to do.' );
    case 'Gaussian'
        logPr.obs = calcMargLogPrData_GaussianNormalInvWishart( Xstats, model.obsModel.params );
    case 'AR-Gaussian'
        logPr.obs = calcMargLogPrData_ARMatrixNormalInvWishart( Xstats, model.obsModel.params );
end


logPr.all = logPr.obs + logPr.z + logPr.F;

end