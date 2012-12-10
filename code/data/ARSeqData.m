% ARSeqData defines a collection of sequential observations
%   that can be modeled jointly by the BP-HMM using Auto-Regressive
%   likelihood models. 
% This object provides a coherent interface to accessing data,
% This includes (1) joint access for computing suff. stats across all seqs
%               (2) individual seq. access
% The particular nature of AR data requires access to both
%   each sequence's raw observations    ( y(t) in E. Fox's thesis)
%        seq(ii) returns a matrix X, which is D x T
%    and aggregated lagged observations (tilde y in E. Fox's thesis)
%        prev(ii) returns matrix XP, which is (D*R) x T
%           where XP(:,tt) is a col vector that stacks the R previous
%           observations, e.g.  [ X(:,tt-1); X(:,tt-2); ... X(:,tt-R) ]
%                  
classdef ARSeqData < SeqData
    properties
        XprevR;
        R;
    end
    methods
        function obj = ARSeqData( varargin )
            obj.R = varargin{1};                
            if length( varargin ) >= 2
                if isobject( varargin{2} )
                    Q = varargin{2};
                else
                    Q = obj@SeqData( varargin(2:end) );
                end
                N = Q.N;
                
                for ii = 1:N
                    if isprop( Q, 'zTrueAll') && ~isempty( Q.zTrueAll )
                        obj = addSeq( obj, Q.seq(ii), Q.name(ii), Q.zTrue(ii) );
                    else
                        obj = addSeq( obj, Q.seq(ii), Q.name(ii) );
                    end
                end
                
            end
            
        end
        
        function retStr = toString(obj )
           retStr = sprintf( 'Emissions Type: %s. Dim: %d. ARorder: %d', obj.getObsType(), obj.D, obj.R ); 
        end
        
        function zTrue = zTrue(obj, ii )
            zTrue = obj.zTrueAll(  obj.aggTs(ii)+1:obj.aggTs(ii+1) );
        end
        
        function Xseq = seq( obj, ii )
            Xseq= obj.Xdata( :, obj.aggTs(ii)+1:obj.aggTs(ii+1) );
        end
        
        function Xprev = prev(obj, ii)
            Xprev=obj.XprevR( :, obj.aggTs(ii)+1:obj.aggTs(ii+1) );
        end
        
        % Add particular sequence's observed data to this collection
        %   Ensure that the data is preprocessed appropriately
        %     so that we obtain the correct AR stats.
        function obj = addSeq( obj, seqData, name, zTrue )
            R = obj.R;
            T=size( seqData, 2 )-R;
            obj.N = obj.N+1;
            obj.D = size( seqData,1);
            obj.Ts(end+1) = T;
            obj.aggTs(end+1) = obj.aggTs(end)+T;
            obj.seqNames{end+1} = name;
            obj.Xdata( :, obj.aggTs(end-1)+1:obj.aggTs(end) ) = seqData(:,R+1:end);
            Xprev = seqData(:,1:end-1);
            XprevR = zeros( obj.D*obj.R, T);
            for tt = 1:T
               XprevR(:,tt) = reshape( Xprev(:,tt+R-1:-1:tt), 1, obj.D*obj.R );
            end
            obj.XprevR(:, obj.aggTs(end-1)+1:obj.aggTs(end) ) = XprevR;
            if exist( 'zTrue','var')
                obj.zTrueAll( obj.aggTs(end-1)+1:obj.aggTs(end) )    = zTrue(R+1:end);
            end
        end
        
    end
end
