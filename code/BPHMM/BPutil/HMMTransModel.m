%HMMTransModel
%  represents HMM transition parameters of many sequences at once
%  provides all necessary methods for sampling, suff. statistics, etc.
classdef HMMTransModel
    properties
        N;
        K;
        seq;
        prior;
        Zstats;
    end
    methods
        function obj = HMMTransModel( N, K)
            obj.N = N;
            obj.K = K;
            obj.seq = repmat( struct('availFeatIDs',[],'eta', [] ), 1, N);
        end
        
        % Set the sticky transition prior used for each sequence
        %   alpha ~ Gam( a_alpha, b_alpha )
        %   kappa ~ Gam( a_kappa, b_kappa )
        %   eta   ~ Gam( alpha + delta*kappa )
        function obj = setPrior( obj, alpha, kappa, a_alpha, b_alpha, a_kappa, b_kappa )
            obj.prior.alpha = alpha;
            obj.prior.kappa = kappa;
            obj.prior.a_alpha = a_alpha;
            obj.prior.b_alpha = b_alpha;
            obj.prior.a_kappa = a_kappa;
            obj.prior.b_kappa = b_kappa;
        end
        
        
        % reallocateFeatIDs( keepFeatIDs )
        %   When reorganizing global features, we need to swap/relabel
        % INPUT: keepFeatIDs
        %           col vector of IDs to keep (need not be in rank order)
        function obj = reallocateFeatIDs( obj, curF, keepFeatIDs )
            obj.K = length( keepFeatIDs );
            for ii = 1:obj.N
                f_ii = false( 1, size(curF,2) );
                availFeatIDs = obj.seq(ii).availFeatIDs;
                f_ii( availFeatIDs ) = 1;
                assert( all( f_ii == curF(ii,:) ), 'Bad feature allocation!' );
                sortAvailIDs = zeros( size( availFeatIDs ) );
                for kk = 1:length(availFeatIDs)
                    sortAvailIDs(kk) = find( keepFeatIDs == availFeatIDs(kk) );
                end
                [newFeatIDs, sortAvailIDs] = sort( sortAvailIDs );
                obj.seq(ii).availFeatIDs = newFeatIDs;
                obj.seq(ii).eta = obj.seq(ii).eta(sortAvailIDs,sortAvailIDs);
            end
        end
        
        
        %setEta()
        %  Deterministically fix eta for particular sequence "ii"
        %    to the given values
        function obj = setEta(obj, ii, f_ii, eta)
           obj.seq(ii).availFeatIDs = find( f_ii );
           obj.seq(ii).eta = eta;
        end
        
        %getSeqEtaWithFeats()
        %  Retrieve the eta parameters for a sequence,
        %    excluding all features not in the given feature vec f_ii
        function seqEta = getSeqEtaWithFeats( obj, ii, f_ii )
           seqEtaIN = obj.seq(ii);
           seqEta = seqEtaIN;
           ksIN = find( f_ii );
           for cc = 1:length( ksIN );
              jjs(cc) = find( seqEtaIN.availFeatIDs == ksIN(cc) );
           end
           seqEta.availFeatIDs = seqEtaIN.availFeatIDs( jjs );
           seqEta.eta = seqEtaIN.eta( jjs, jjs );
        end
        
        %getAllEta_PriorMean()
        %  Set the eta for every sequence to its prior mean
        function obj = getAllEta_PriorMean( obj, seqIDs, F, propFeatIDs )
          obj.K = size(F,2);
          propON = false(1,obj.K);
          propON(propFeatIDs) = 1;
          for ii = seqIDs
             kIDs = find( F(ii,:) | propON );
             obj.seq(ii).availFeatIDs = kIDs;
             Kii = length( obj.seq(ii).availFeatIDs );
             obj.seq(ii).eta = obj.prior.alpha*ones(Kii,Kii) + obj.prior.kappa*eye(Kii);
          end
        end
        
        % Deterministically set eta of all sequences to given values
        % INPUT:
        %   F : N x K feature matrix
        %   EtaAll : N x1 cell array, where EtaAll{ii} is K x K matrix
        function obj = updateAllEta( obj, F, EtaAll, seqIDs )
           if ~exist( 'seqIDs', 'var' )
               seqIDs = 1:size(F,1);
           end
           for ii = seqIDs
              kIDs = find( F(ii,:) );
              obj.seq(ii).availFeatIDs = kIDs;
              obj.seq(ii).eta = EtaAll{ii}( kIDs, kIDs ); 
           end
        end
        
        function obj = updateAllZSuffStats( obj, F, stateSeq )
            obj.Zstats = struct('Nz', [], 'Nz_init', [], 'availFeatIDs',  []);
            if isempty( stateSeq ) % just set to all zeros
                for ii = 1:size(F,1)
                    obj.Zstats(ii) = obj.getZSuffStats( F(ii,:), [] );
                end
            else
                for ii = 1:size(F,1)
                    obj.Zstats(ii) = obj.getZSuffStats( F(ii,:), stateSeq(ii).z  );
                end
            end
        end
        
        function Zstats = getZSuffStats( obj, f_ii, z_ii )
           availFeatIDs = find( f_ii );
           Kii = length( availFeatIDs );
           Zstats.Nz      = zeros( Kii, Kii );
           Zstats.Nz_init = zeros( 1, Kii );
           Zstats.availFeatIDs = availFeatIDs;
           if ~isempty( z_ii )
               zprev = z_ii(1:end-1);
               zcur  = z_ii(2:end);
               Zstats.Nz_init = (  z_ii(1)==availFeatIDs );
               for jj = 1:Kii
                   Zstats.Nz(jj,:) = histc( zcur( zprev==availFeatIDs(jj)  ), availFeatIDs );
               end
           end
        end
        
        function eta = sampleEta( obj, f_ii, z_ii )
           Zstats = obj.getZSuffStats( f_ii, z_ii );
           eta = obj.sampleEta_FromZStats( Zstats );
        end
        
        function obj = sampleAllEta( obj, F, stateSeq )
            if ~exist( 'stateSeq', 'var' )
                stateSeq = [];
            end
           obj.K = size(F,2);
           obj = obj.updateAllZSuffStats( F, stateSeq );
           for ii = 1:size(F,1)
              obj.seq(ii).availFeatIDs = obj.Zstats(ii).availFeatIDs;
              obj.seq(ii).eta = obj.sampleEta_FromZStats( obj.Zstats(ii) ); 
           end
        end
        
        function pi = pi( obj, ii )
           pi = obj.seq(ii).eta;
           pi = bsxfun( @rdivide, pi, sum(pi,2) );
        end
        
        function pi0 = pi_init( obj, ii )
            K = length(obj.seq(ii).availFeatIDs );
            pi0 = ones( 1,  K )./K;
        end
        
        function eta = sampleEta_FromZStats(obj, Zstats )
            alpha = obj.prior.alpha;
            kappa = obj.prior.kappa;
            Kii   = size( Zstats.Nz, 1 );
            % Draw Normalized Probabilities from Dirichlet Posterior
            %   q ~ Dir( N_k + a + kappa*delta(j,k) )
            Pi_z  = randgamma(  Zstats.Nz+ alpha*ones(Kii,Kii) + kappa*eye( Kii,Kii ) );
            Pi_z  = bsxfun( @rdivide, Pi_z, sum( Pi_z,2) );
            % Draw a scale factor for each row of Pi_z
            %   proportional to *sum* of prior parameters
            EtaSum = randgamma( (kappa + Kii*alpha)*ones(Kii,1) );
            % Combine Dir draws with scale factor to get Gamma draws
            % via the transformation:
            %    eta_k = q_k * EtaSum  where   sum( eta_k ) = EtaSum
            eta = bsxfun( @times, Pi_z, EtaSum );
        end
        
        function propEta = sampleEtaProposal_Shared(obj, ii)
            alpha = obj.prior.alpha;
            kappa = obj.prior.kappa;
            
            propEta = randgamma( alpha*ones(obj.K,obj.K) + kappa*eye(obj.K,obj.K) );
            ks = obj.seq(ii).availFeatIDs;
            propEta(ks,ks) = obj.seq(ii).eta;
            
        end
        
        
        function Eta = sampleEtaProposal_UniqueBirth(obj, ii )
            alpha = obj.prior.alpha;
            kappa = obj.prior.kappa;
            
            Kii = length( obj.seq(ii).availFeatIDs );
            newEtaCol = randgamma( [alpha*ones(Kii, 1); alpha+kappa] );
            newEtaRow = randgamma( alpha*ones(1, Kii ) );
            Eta = obj.seq(ii).eta;
            Eta(1:Kii+1, Kii+1) = newEtaCol;
            Eta(Kii+1, 1:Kii) = newEtaRow;
        end
        
        function logPr = calcMargPrStateSeq( obj, F, stateSeq, objIDs )
            if ~exist('objIDs','var')
                objIDs = 1:size(F,1);
            end
            obj = updateAllZSuffStats( obj, F, stateSeq );
            logPrZ = zeros(1, obj.N);
            for ii = objIDs
                logPrZ(ii) = calcMargLogPrData_MultinomialDirichletSticky( obj.Zstats(ii).Nz, obj.prior.alpha, obj.prior.kappa );
            end
            logPr = sum( logPrZ );
        end
        
        
    end
end
