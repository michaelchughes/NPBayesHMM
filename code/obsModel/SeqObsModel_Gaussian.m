classdef SeqObsModel_Gaussian < SeqObsModel
    % SeqObsModel : defines Gaussian emission model for observed data in HMM
    %   Allows for multivariate observations of any dimension.
    %   Assumes conjugate prior on mean and covariance (Normal-Inverse Wishart).
    %   Extends the base object "SeqObsModel", filling in Gaussian details
    methods
        % ====================================================== CONSTRUCTOR
        function obj = SeqObsModel_Gaussian( K, D)
            obj.K = K;
            obj.D = D;
            obj.theta = repmat( struct('mu',zeros(D,1),'invSigma', zeros(D,D) ), 1, K );
        end
        
        % setPrior() sets prior on emission parameters to specified values.
        % USAGE:
        %   To set to default values given data object and model params
        %     ObsM = ObsM.setPrior( data, modelParamStruct )
        %     where modelParamStruct has attributes:
        %        degFree : # deg of freedom
        %        precMu  : precision parameter for mean
        %        Scoef   : scalar multiplier for Scale Matrix
        %        doEmpCovScalePrior: boolean.
        %               true  = set ScaleMat to Scoef*emp covariance
        %               false = set ScaleMat to Scoef*eye(D)
        %  To initialize directly to given parameters
        %      ObsM = ObsM.setPrior( mu, precMu, degFree, ScaleMat )
        function obj = setPrior( obj, varargin )
            obj.priorDef = struct();
            obj.priorDef.doEmpCovScalePrior =0;
            obj.priorDef.Scoef = 0;
            if strfind( class(varargin{1} ), 'SeqData' )
                data = varargin{1};
                obsM = varargin{2};
                obj.prior.mu = zeros( data.D, 1 );
                obj.prior.precMu = obsM.precMu;
                obj.prior.degFree = max( obsM.degFree, data.D+2 );
                obj.priorDef.doEmpCovScalePrior = obsM.doEmpCovScalePrior;
                obj.priorDef.Scoef = obsM.Scoef;
                if obsM.doEmpCovScalePrior
                    obj.prior.ScaleMat = obsM.Scoef * cov( data.Xdata', 1);
                else
                    obj.prior.ScaleMat = obsM.Scoef * eye( data.D );
                end
            elseif length( varargin ) == 1
                PStruct = varargin{1};
                obj.prior.mu = PStruct.Mu(:);
                obj.prior.precMu = PStruct.precMu;
                obj.prior.degFree = PStruct.degFree;
                obj.prior.ScaleMat = PStruct.ScaleMat;
            else
                obj.prior.mu = varargin{1}(:);
                obj.prior.precMu = varargin{2};
                obj.prior.degFree = varargin{3};
                obj.prior.ScaleMat = varargin{4};
            end            
        end
        
        % ==================================================== GET Xstats
        function Xstats = getXSuffStats( obj, Xkk )
            Xstats.nObs = size(Xkk,2);
            if Xstats.nObs > 0
                Xstats.Xsum = sum( Xkk, 2);
                Xstats.XXsum = Xkk*Xkk';
            else
                Xstats.Xsum = [];
                Xstats.XXsum = [];
            end
        end
        
        function obj = incXStats( obj, ii, data, stateSeq, featIDs )
            for kk = featIDs
                Xkk = data.seq(ii);
                Xkk = Xkk(:, stateSeq(ii).z == kk );
                nNew = size(Xkk,2);
                if length(obj.Xstats)<kk || obj.Xstats(kk).nObs == 0
                    obj.Xstats(kk).nObs = nNew;
                    obj.Xstats(kk).Xsum = sum(Xkk,2);
                    obj.Xstats(kk).XXsum = Xkk*Xkk';
                else
                    obj.Xstats(kk).nObs = obj.Xstats(kk).nObs + nNew;
                    obj.Xstats(kk).Xsum = obj.Xstats(kk).Xsum + sum(Xkk,2);
                    obj.Xstats(kk).XXsum = obj.Xstats(kk).XXsum + Xkk*Xkk';
                end
            end
        end
        
        function obj = decXStats( obj, ii, data, stateSeq, featIDs )
            for kk = featIDs
                Xkk = data.seq(ii);
                Xkk = Xkk(:, stateSeq(ii).z == kk );
                nNew = size(Xkk,2);
                if obj.Xstats(kk).nObs-nNew <= 0
                    obj.Xstats(kk).nObs=0;
                    obj.Xstats(kk).Xsum=[];
                    obj.Xstats(kk).XXsum=[];
                else
                    obj.Xstats(kk).nObs = obj.Xstats(kk).nObs - nNew;
                    obj.Xstats(kk).Xsum = obj.Xstats(kk).Xsum - sum(Xkk,2);
                    obj.Xstats(kk).XXsum = obj.Xstats(kk).XXsum - Xkk*Xkk';
                end
            end
        end
        
        % ==================================================== GET params
        function retStr = getParamDescr( obj, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            retStr = 'Norm-InvWish. Mu0=0, dFree=%d,';
            if obj.priorDef.doEmpCovScalePrior
                retStr = [retStr ' S0=%.2f*EmpCov.'];
            else                
                retStr = [retStr ' S0=%.2f*eye.'];
            end
            retStr = sprintf( retStr, obj.prior.degFree, obj.priorDef.Scoef );
        end
        
        function PP = getPosteriorParams( obj, Xstats )
            if Xstats.nObs > 0
                N = Xstats.nObs;
                degFreeN = N + obj.prior.degFree;
                precMuN  = N + obj.prior.precMu;
                Xdiff = Xstats.Xsum/N - obj.prior.mu;
                XdiffMat = ( obj.prior.precMu * N  )/( obj.prior.precMu + N ) .* ( Xdiff*Xdiff' );
                CovMat = (Xstats.XXsum - Xstats.Xsum*Xstats.Xsum'/N);
                ScaleMatN = obj.prior.ScaleMat + CovMat + XdiffMat;
                MuN = (Xstats.Xsum + obj.prior.mu*obj.prior.precMu)/ ( precMuN );
                
                PP.mu = MuN;
                PP.precMu = precMuN;
                PP.degFree = degFreeN;
                PP.ScaleMat = ScaleMatN;
            else
                PP = obj.prior;
            end
        end
        
        
        
        % ==================================================== GET theta mean
        
        % Retrieve the expected value of theta under the given params
        % In Gaussian case,
        %   mu = PP.mu
        %   Sigma = PP.ScaleMat / (degFree-D-1 )
        function [theta] = getTheta_Mean( obj, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            theta.mu  = PP.mu;
            theta.invSigma = (PP.degFree-obj.D-1)*( PP.ScaleMat \ eye(obj.D) );
        end
        
        % ==================================================== GET theta samps
        function [theta, PP] = sampleTheta_FromParams( obj, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            [~, sqrtInvSigma] = randiwishart( PP.ScaleMat, PP.degFree );
            theta.mu = ( sqrt( PP.precMu ) .* sqrtInvSigma ) \randn(obj.D,1)  + PP.mu;
            theta.invSigma = sqrtInvSigma'*sqrtInvSigma;
        end
        
        
        function [theta, PN] = sampleTheta( obj, Xkk )
            if ~exist( 'Xkk', 'var' )
                Xkk = [];
            end
            Xstats = obj.getXSuffStats( Xkk );
            PN = obj.getPosteriorParams( Xstats );
            [theta] = obj.sampleTheta_FromParams( PN );
        end
        
        % ==================================================== GET log probs
        
        % Calc probability of a realization of a emission parameter "theta"
        %   given the provided prior parameter struct "PP"
        function logPr = calcLogPrTheta( obj, theta, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            logPr = calcLogPrNormalInvWishart( theta.mu, theta.invSigma, PP );
        end
        
        % Calculate soft evidence for particular sequence "ii"
        %   given all available emission parameters (indicated by "kIDs")
        %OUTPUT
        %  logSoftEv : K x T matrix, where entry k,t gives log probability
        %                of the t-th observation under theta(kk)
        function logSoftEv = calcLogSoftEv(obj, ii, data, kIDs )
            if ~exist( 'kIDs', 'var' )
                kIDs = 1:obj.K;
            end
            Xseq = data.seq(ii);
            T = size( Xseq,2);
            logSoftEv = -inf( obj.K, T );
            for kk = kIDs
                cholInvSigma = chol( obj.theta(kk).invSigma );
                logDetInvSigma = 2*sum( log( diag( cholInvSigma) ) );
                XdiffMu = bsxfun(@minus, Xseq, obj.theta(kk).mu );
                U = XdiffMu'*cholInvSigma';
                logSoftEv(kk,:) = 0.5*logDetInvSigma - 0.5*sum( U.^2,2);
            end
            logSoftEv = logSoftEv - 0.5*obj.D*log(2*pi);
        end
        
        % Calculate marginal probability of observed data
        %  given the state sequence assignments of every observation
        % USAGE:
        %  logPr = ThetaM.calcMargPrData()  
        %    calculate using *stored* sufficient statistics in this object 
        %  logPr = ThetaM.calcMargPrData( data, stateSeq )
        %    calculate using freshly computed suff stats
        function logPr = calcMargPrData(obj, data, stateSeq, ks )
            if ~exist( 'ks', 'var' )
                ks = 1:obj.K;
            end
            if exist( 'data', 'var' ) && ~isempty( data )
                obj = obj.updateAllXSuffStats(  horzcat(stateSeq(:).z), data, ks);
            end
            LOG_PI =  1.144729885849400;
            P0 = obj.prior;
            D = obj.D;
            logPr = zeros(1, obj.K);
            for kk = ks
                Nkk = obj.Xstats(kk).nObs;
                if Nkk > 0
                    PN = obj.getPosteriorParams( obj.Xstats(kk) );
                    
                    logPr(kk) = ...
                        logMvGamma(0.5*PN.degFree,D) - logMvGamma(0.5*P0.degFree,D) ...
                        + 0.5*P0.degFree*log( det( P0.ScaleMat ) )   ...
                        - 0.5*PN.degFree*log( det( PN.ScaleMat ) ) ...
                        + 0.5*D*log( P0.precMu / PN.precMu );
                end
            end
            N = sum( [obj.Xstats(ks).nObs] );
            logPr = sum(logPr) - 0.5*N*D*LOG_PI;
        end
        
    end
end
