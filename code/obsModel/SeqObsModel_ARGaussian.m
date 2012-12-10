classdef SeqObsModel_ARGaussian < SeqObsModel
    % SeqObsModel : defines autoregressive emission model for data in HMM
    %   Allows for multivariate observations of any dimension (param 'R')
    %   Allows for any # of autoregressions, eg AR(1), AR(2), ... AR(R)
    %   Assumes conjugate prior on mean and covariance (Normal-Inverse Wishart).
    %   Extends the base object "SeqObsModel", filling in Gaussian details
    properties
        R;
    end
    
    methods
        % ====================================================== CONSTRUCTOR
        function obj = SeqObsModel_ARGaussian( K, D, R)
            obj.K = K;
            obj.D = D;
            obj.R = R;
            obj.theta = repmat( struct('A',zeros(D,D*R),'invSigma', zeros(D,D) ), 1, K );
        end
        
        % setPrior() sets prior on emission parameters to specified values.
        % USAGE:
        %   To set to default values given data object and model params 
        %     ObsM = ObsM.setPrior( data, modelParamStruct )
        %     where modelParamStruct has attributes:
        %        degFree : # deg of freedom
        %        Scoef   : scalar multiplier for Scale Matrix
        %        doEmpCovScalePrior: boolean. 
        %               true  = set ScaleMat to Scoef*emp covariance
        %               false = set ScaleMat to Scoef*eye(D)
        %  To initialize directly to given parameters
        %      ObsM = ObsM.setPrior( MeanMat, invAScaleMat, degFree, ScaleMat )
        function obj = setPrior( obj, varargin )
            obj.priorDef.doEmpCov = 0;
            obj.priorDef.doEmpCovFirstDiff = 0;
            obj.priorDef.Scoef = 0;
            if strfind( class(varargin{1} ), 'SeqData' )
                data = varargin{1};
                obsM = varargin{2};                
                obj.prior.MeanMat = zeros( data.D, data.D*obj.R );
                obj.prior.invAScaleMat = obsM.Scoef*eye( data.D*obj.R, data.D*obj.R);
                obj.prior.degFree = max( obsM.degFree, data.D+2 );
                obj.priorDef.doEmpCov = obsM.doEmpCov;
                obj.priorDef.doEmpCovFirstDiff = obsM.doEmpCovFirstDiff;
                obj.priorDef.Scoef  = obsM.Scoef;                
                if obsM.doEmpCov
                    obj.prior.ScaleMat = obsM.Scoef * cov( data.Xdata', 1);
                elseif obsM.doEmpCovFirstDiff
                    obj.prior.ScaleMat = obsM.Scoef * cov( diff(data.Xdata'), 1); 
                else
                    obj.prior.ScaleMat = obsM.Scoef * eye( data.D );
                end
            elseif length( varargin ) == 1
                PStruct = varargin{1};
                obj.prior.MeanMat = PStruct.MeanMat;
                obj.prior.invAScaleMat = PStruct.invAScaleMat;
                obj.prior.degFree = PStruct.degFree;
                obj.prior.ScaleMat = PStruct.ScaleMat;
            else
                obj.prior.MeanMat = varargin{1}(:);
                obj.prior.invAScaleMat = varargin{2};
                obj.prior.degFree = varargin{3};
                obj.prior.ScaleMat = varargin{4};
            end
        end
        
        % ==================================================== GET Xstats
        function Xstats = getXSuffStats( obj, Xkk, XkkPrev )
            Xstats.nObs = size(Xkk,2);
            if Xstats.nObs > 0
                Xstats.XX = Xkk*Xkk';
                Xstats.XY = Xkk*XkkPrev';
                Xstats.YY = XkkPrev*XkkPrev';                
            else
                Xstats.XX = [];
                Xstats.XY = [];
                Xstats.YY = [];
            end
        end
        
          function obj = incXStats( obj, ii, data, stateSeq, featIDs )
            for kk = featIDs                
                Xkk = data.seq(ii);
                Xprev = data.prev(ii);
                
                kkINDS = stateSeq(ii).z==kk;
                Xkk = Xkk(:, kkINDS );
                Xprev = Xprev(:, kkINDS );
                
                nNew = size(Xkk,2);
                XX = Xkk*Xkk';
                XY = Xkk*Xprev';
                YY = Xprev*Xprev';
                if length(obj.Xstats)<kk ||obj.Xstats(kk).nObs == 0
                    obj.Xstats(kk).nObs = nNew;
                    obj.Xstats(kk).XX = XX;
                    obj.Xstats(kk).XY = XY;
                    obj.Xstats(kk).YY = YY;
                else
                    obj.Xstats(kk).nObs = obj.Xstats(kk).nObs + nNew;
                    obj.Xstats(kk).XX = obj.Xstats(kk).XX + XX;
                    obj.Xstats(kk).XY = obj.Xstats(kk).XY + XY;
                    obj.Xstats(kk).YY = obj.Xstats(kk).YY + YY;
                end
            end
        end
        
        function obj = decXStats( obj, ii, data, stateSeq, featIDs )
            for kk = featIDs
                Xkk = data.seq(ii);
                Xprev = data.prev(ii);
                
                kkINDS = stateSeq(ii).z==kk;
                Xkk = Xkk(:, kkINDS );
                Xprev = Xprev(:, kkINDS );
                
                nNew = size(Xkk,2);
                XX = Xkk*Xkk';
                XY = Xkk*Xprev';
                YY = Xprev*Xprev';
                
                if obj.Xstats(kk).nObs-nNew <= 0
                    obj.Xstats(kk).nObs=0;
                    obj.Xstats(kk).XX=[];
                    obj.Xstats(kk).XY=[];
                    obj.Xstats(kk).YY=[];
                else
                    obj.Xstats(kk).nObs = obj.Xstats(kk).nObs - nNew;                    
                    obj.Xstats(kk).XX = obj.Xstats(kk).XX - XX;
                    obj.Xstats(kk).XY = obj.Xstats(kk).XY - XY;
                    obj.Xstats(kk).YY = obj.Xstats(kk).YY - YY;
                end
            end
        end

        % ==================================================== GET params
        function retStr = getParamDescr( obj, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            retStr = 'MNIW. dFree=%d,';
            if obj.priorDef.doEmpCov
                retStr = [retStr ' S0=%.2f*EmpCov'];
            elseif obj.priorDef.doEmpCovFirstDiff                
                retStr = [retStr ' S0=%.2f*FirstDiffEmpCov'];
            else
                retStr = [retStr ' S0=%.2f*eye'];            
            end
            retStr = sprintf( retStr, obj.prior.degFree, obj.priorDef.Scoef );
        end

        function PP = getPosteriorParams( obj, Xstats )
            if Xstats.nObs > 0
                degFreeN = Xstats.nObs + obj.prior.degFree;
                
                MK  = obj.prior.MeanMat*obj.prior.invAScaleMat;
                MKM = MK*obj.prior.MeanMat';
                Sxx = Xstats.XX + MKM;
                Syy = Xstats.YY + obj.prior.invAScaleMat;
                Sxy = Xstats.XY + MK;
                
                invAScaleMatN = Syy;
                MeanMatN  = Sxy / Syy; % Sxy * inv(Syy)
                ScaleMatN = obj.prior.ScaleMat + Sxx - (MeanMatN)*Sxy';
                            
                PP.MeanMat = MeanMatN;
                PP.invAScaleMat = invAScaleMatN;
                PP.degFree = degFreeN;
                PP.ScaleMat = ScaleMatN;
            else
                PP = obj.prior;
            end
        end
        
         % ==================================================== GET theta mean
        function [theta] = getTheta_Mean( obj, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            theta.A  = PP.MeanMat;
            theta.invSigma = (PP.degFree-obj.D-1)*( PP.ScaleMat \ eye(obj.D) );
        end

        % ==================================================== GET theta samps
        function [theta, PP] = sampleTheta_FromParams( obj, PP )
           if ~exist('PP','var') 
               PP = obj.prior;
           end
            [sqrtSigma, sqrtInvSigma] = randiwishart( PP.ScaleMat, PP.degFree );
            theta.A = sampleFromMatrixNormal( PP.MeanMat, sqrtSigma, chol( inv(PP.invAScaleMat) ) );
            theta.invSigma = sqrtInvSigma'*sqrtInvSigma;
        end
                
        
        function [theta, PN] = sampleTheta( obj, Xkk, Xprev )
            if ~exist( 'Xkk', 'var' )
                Xkk = [];
                Xprev = [];
            end
            Xstats = obj.getXSuffStats( Xkk, Xprev );
            PN = obj.getPosteriorParams( Xstats );
            [theta] = obj.sampleTheta_FromParams( PN );
        end
        
        % ==================================================== GET log probs
        function logPr = calcLogPrTheta( obj, theta, PP )
            if ~exist('PP','var')
                PP = obj.prior;
            end
            logPr = calcLogPrMatrixNormalInvWishart( theta.A, theta.invSigma, PP );
        end

        function logSoftEv = calcLogSoftEv(obj, ii, data, kIDs )
            if ~exist( 'kIDs', 'var' )
                kIDs = 1:obj.K;
            end
            Xseq = data.seq(ii);
            Xprev = data.prev(ii);
            
            T = size( Xseq,2);
            logSoftEv = -inf( obj.K, T );
            for kk = kIDs
                cholInvSigma  = chol( obj.theta(kk).invSigma );
                logDetInvSigma = 2*sum( log( diag( cholInvSigma ) ) );

                XdiffMu = Xseq - obj.theta(kk).A*Xprev;
                U = XdiffMu'*cholInvSigma';
                
                logSoftEv(kk,:) = 0.5*logDetInvSigma - 0.5*sum( U.^2,2);
            end
            logSoftEv = logSoftEv - 0.5*obj.D*log(2*pi);
        end
        
        
        function logPr = calcMargPrData(obj, data, stateSeq, ks)
            if ~exist( 'ks', 'var' )
                ks = 1:obj.K;
            end
            if exist( 'data', 'var' ) && ~isempty(data)
                obj = obj.updateAllXSuffStats( horzcat(stateSeq(:).z), data, ks );
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
                      + logMvGamma( 0.5*PN.degFree, D ) ...
                      - logMvGamma( 0.5*P0.degFree, D ) ...
                      + 0.5*P0.degFree*log( det( P0.ScaleMat ) )   ...
                      - 0.5*PN.degFree*log( det( PN.ScaleMat ) ) ...
                      + 0.5*D*log( det( P0.invAScaleMat ) ) ...
                      - 0.5*D*log( det( PN.invAScaleMat ) );
                end
            end
            N = sum( [obj.Xstats(ks).nObs] );
            logPr = sum(logPr) - 0.5*N*D*LOG_PI;
        end
        
    end
end
