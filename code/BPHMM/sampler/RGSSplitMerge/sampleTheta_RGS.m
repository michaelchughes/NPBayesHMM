function [theta, logQ_RGS] = sampleTheta_RGS( F, stateSeq, theta, data_struct, model, selectFeatIDs, TargetPsi )
% Sample new value for emission parameters Theta under RESTRICTED setting
%     at *select* features only
% Output theta struct has same params as input struct
%    except for features in vector "selectFeatIDs", 
%          which are brand new posterior draws
% TO DO: Speed up by only getting sufficient stats for active objs/feats
%activeObjIDs = find( F(:,selectFeatIDs) );

logQ_RGS = 0;

[Xstats] = getXSuffStats( F, stateSeq, data_struct, model );

switch model.obsModel.type
    % ---------------------------------------------------------- Multinomial     
    case 'Multinomial'
       [K V] = size( Xstats.Nkv );
       lam_vec = model.obsModel.params.lambda;

       Kselect = length( selectFeatIDs );
       for jj = 1:Kselect
          kz = selectFeatIDs( jj );
          Nv =[ lam_vec+Xstats.Nkv(kz,:)'  ] ;
          
          if exist( 'TargetPsi','var' ) && ~isempty( TargetPsi )
            matchFeatID = TargetPsi.activeFeatIDs( jj );
            theta.logp(kz,:) = TargetPsi.theta.logp( matchFeatID,: );
          else
            theta.logp(kz,:) = log( randdirichlet(  Nv )' );
          end
          
          if nargout > 1
              logQ_RGS = logQ_RGS + calcLogPrDirichlet( theta.logp(kz,:), Nv', 1 );
          end
       end
       
    % --------------------------------------------------------  Bernoulli   
    case 'Bernoulli'
        error('TO DO');
       [K D] = size( Xstats.nON );
       a_p = model.obsModel.params.a_p;
       b_p = model.obsModel.params.b_p;
       
       P = zeros( K, D );
       for kz = 1:K
           A  = Xstats.nON(kz,:) + a_p;
           B  = Xstats.nOFF(kz,:) + b_p;
           P(kz,:) = betarnd( A,  B );
       end
       theta.p = P;
       
    % --------------------------------------------------------  Gaussian
    case 'Gaussian'
        [K D] = size( Xstats.Mu );
        
        Kselect = length( selectFeatIDs );
           
        PP = model.obsModel.params;
        degFree = PP.degFree + Xstats.nObs;
        precMu  =  PP.precMu + Xstats.nObs;
        
        switch model.obsModel.priorType
            case 'NormalInvWishart'
                ScaleMatrix = zeros(D, D, K);
                for jj = 1:Kselect
                    kz = selectFeatIDs( jj );
                    Sjj = Xstats.Mu(kz,:) - PP.Mu;
                    A = ( PP.precMu * Xstats.nObs(kz)  )/( PP.precMu + Xstats.nObs(kz)  ) .* ( Sjj' * Sjj );
                    ScaleMatrix(:,:,kz)  = PP.ScaleMatrix + Xstats.Sigma(:,:,kz) + A;
                end
                MuN = bsxfun( @plus, bsxfun(@times, Xstats.Mu, Xstats.nObs), PP.Mu * PP.precMu );
                MuN = bsxfun( @rdivide, MuN, Xstats.nObs+PP.precMu );
                
                for jj = 1:Kselect
                    kz = selectFeatIDs( jj );
                    [~, sqrtInvSigma] = randiwishart( ScaleMatrix(:,:,kz), degFree(kz) );
                    invSjj = sqrtInvSigma'*sqrtInvSigma;
                    theta.invSigma(:,:,kz) = invSjj;
                    theta.Mu(kz,:)     = ( sqrt(precMu(kz)) .* sqrtInvSigma ) \randn(D,1)  + MuN(kz,:)';
                end                
        end
        
    % --------------------------------------------------------  AR   
    case 'AR-Gaussian' 
        R = model.obsModel.ARorder;                
        D = size( data_struct(1).obs, 2 );                
        K = size( F, 2 );
        Kselect = length( selectFeatIDs );

        PP = model.obsModel.params;
        CC = PP.invColScaleMatrix;
        
        for jj = 1:Kselect
            kz = selectFeatIDs( jj );
            %  XX : data.obs
            %  YY : data.obsPrevR
            XX = Xstats.XX(:,:,kz);
            XY = Xstats.XY(:,:,kz);
            YY = Xstats.YY(:,:,kz);
            
            % ------------------------------------  draw Sigma                        
            degFree_jj = PP.degFree + Xstats.nObs(kz);
            S = XX - ( XY/(YY + CC) )*XY';
            %S = (S + S')/2; % enforce symmetry numerically
            SM_jj = PP.ScaleMatrix + S;

            [sqrtSigma, sqrtInvSigma] = randiwishart( SM_jj, degFree_jj );
            theta.invSigma(:,:,kz) = sqrtInvSigma'*sqrtInvSigma;
            
            %Sigma_jj = sqrtSigma'*sqrtSigma;
            %theta.Sigma(:,:,jj) = Sigma_jj;
            
            % ------------------------------------  draw A | Sigma
            sqrtInvCC_N = chol( inv( YY + CC ) );
            M_N         = XY / ( YY + CC );  % XY * inv( YY+CC) 
            theta.A(:,:,kz) = sampleFromMatrixNormal( M_N, sqrtSigma, sqrtInvCC_N);
        end

end

end % main function