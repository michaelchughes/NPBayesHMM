classdef SeqObsModel
    % SeqObsModel : defines probabilistic model for observed data in HMM
    % This generic class provides common operations for posterior inference
    %  further specified by descendants, such as "SeqObsModel_Gaussian"
    %  which implement specifics of a particular emission distrib. (Gaussian).
    properties
        K;
        D;
        theta;
        prior;
        priorDef;
        Xstats;
        name;
    end
    
    % --------------------------------------------------- ABSTRACT
    % Any subclass needs to implement all these methods
    %function obj = setPrior( obj, varargin );        
    %function Xstats = getXSuffStats(obj, Xkk );
        
    %function str = getParamDescr( obj, PP ); end
        
    %function PN  = getPosteriorParams(obj, Xstats ); end
       
    %function [theta,PP] = sampleTheta_FromParams(obj,PP); end
    %function [theta,PN] = sampleTheta(obj,Xkk); end
        
    %function logPr = calcLogPrTheta(obj, theta,PP); end
    %function logSoftEv = calcLogSoftEv(obj, Xseq, kIDs ); end
    %function logPr = calcMargPrData(obj, data, stateSeq ); end
   
    methods
        % ===================================================== SET methods
        % These alter the stored object, and must be called with syntax
        % MyObj = MyObj.genericSetMethod( params );


        function obj = insertTheta(obj, theta, kk )
            if ~exist( 'kk', 'var' )
                obj.K = obj.K+1;
                kk = obj.K;
            end
            obj.theta(kk) = theta;
        end
        
        function obj = sampleAllTheta(obj, data, stateSeq )
            Zall = horzcat( stateSeq(:).z );
            obj = obj.updateAllXSuffStats( Zall, data );
            for kk = 1:obj.K
               PN = obj.getPosteriorParams( obj.Xstats(kk) );
               obj.theta(kk) = obj.sampleTheta_FromParams( PN ); 
            end
        end
        
        function obj = getAllTheta_PosteriorMean( obj, data, stateSeq, featIDs, zeroIDs )
            if ~exist('featIDs','var')
                featIDs = []; zeroIDs=[];
            end
            Zall = horzcat( stateSeq(:).z );
            obj = obj.updateAllXSuffStats( Zall, data, featIDs, zeroIDs );
            for kk = 1:obj.K
               PN = obj.getPosteriorParams( obj.Xstats(kk) );
               obj.theta(kk) = obj.getTheta_Mean( PN ); 
            end
        end
        
        function obj = updateAllXSuffStats(obj, Zall, data, featIDs, zeroIDs )
            if ~exist( 'featIDs', 'var' ) || isempty(featIDs)
                featIDs = 1:max(Zall);
                obj.Xstats = repmat( obj.getXSuffStats([]), obj.K, 1 ); 
            else
                obj.K = max( [max(Zall), obj.K, max(featIDs)] );
                if isempty( obj.Xstats )
                    obj.Xstats = repmat( obj.getXSuffStats([]), obj.K, 1 ); 
                end
            end
            if exist( 'zeroIDs', 'var' ) && ~isempty(featIDs)
               for zz = zeroIDs
                 obj.Xstats(zz).nObs = 0; 
               end
            end
            if strcmp( class(data), 'SeqData' )
               for kk = featIDs
                    kkINDS = Zall ==kk;
                    Xkk = data.Xdata(:,kkINDS);
                    obj.Xstats(kk) = obj.getXSuffStats( Xkk );
                end
            elseif strcmp( class(data), 'ARSeqData' )
                for kk = featIDs
                    kkINDS = Zall ==kk;
                    Xkk = data.Xdata(:,kkINDS);
                    Xprev = data.XprevR(:,kkINDS);
                    
                    obj.Xstats(kk) = obj.getXSuffStats( Xkk, Xprev);
                end
            end
        end
        
        function obj = reallocateFeatIDs( obj, keepIDs )
            obj.theta = obj.theta( keepIDs );
            obj.K = length( keepIDs );
            if max(keepIDs) <= length(obj.Xstats)
                obj.Xstats = obj.Xstats(keepIDs);
            else
                obj.Xstats = [];
            end
        end


        % ====================================================== GET methods
        
        % calcMargLikRatio_MergeFeats
        %  logPrDiff = logPrProp - logPrCur, where propose = merge feats ki,kj
        %                                          current = keep ki,kj separate
        function [logPrDiff] = calcMargLikRatio_MergeFeats(obj, data, stateSeq, ki, kj )
            logPrCur = obj.calcMargPrData( [], [], [ki, kj] );
            propStateSeq = stateSeq;
            for ii = 1:length(propStateSeq)
                propStateSeq(ii).z( stateSeq(ii).z==kj ) = ki;
            end
            logPrProp = obj.calcMargPrData( data, propStateSeq, ki );
            logPrDiff = logPrProp - logPrCur;
        end
        
        function [theta,PP] = sampleThetaProposal_BirthPrior( obj )
          [theta,PP] = obj.sampleTheta(); 
        end
        
        function [theta, PPmix, choice] = sampleThetaProposal_BirthDataDriven( obj, ii, data, winIDs )
            if strcmp( class(data), 'ARSeqData' )
                X = data.seq(ii);
                Xprev = data.prev(ii);
                PN = obj.getPosteriorParams( obj.getXSuffStats( X(:,winIDs), Xprev(:,winIDs) ) );
            else
                X = data.seq(ii);
                PN = obj.getPosteriorParams( obj.getXSuffStats( X(:,winIDs) ) );
            end
            if rand < 0.5
              choice = 1;
              theta = obj.sampleTheta_FromParams(); % Prior
            else
              choice = 2;
              theta = obj.sampleTheta_FromParams( PN ); % Posterior
            end
            PPmix(1) = obj.prior;
            PPmix(2) = PN;
        end
        
        
        function logPr = calcLogPrTheta_MixWithPrior( obj, theta, PPmix )
            C = length(PPmix );
            logpi = log(1./C);
            logQ=zeros(C,1);
            for aa = 1:C
                logQ(aa) = logpi + obj.calcLogPrTheta( theta, PPmix(aa) );
            end
            logPr = logsumexp( logQ );
        end
        
    end   
end
