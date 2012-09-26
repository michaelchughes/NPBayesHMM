% SeqData defines a collection of sequential observations
%   that can be modeled jointly by the BP-HMM
% This object provides a coherent interface to accessing data,
%   including (1) joint access for computing suff. stats across all seqs
%             (2) individual seq. access via the .seq( ii ) method
% Possible subclass: SparseSeqData,  ARSeqData (for VAR models)
classdef SeqData
    properties
        Xdata;
        N;
        D;
        Ts;
        aggTs;
        zTrueAll;
        seqNames;
    end
    methods

        function obj = SeqData( varargin )
            if(nargin > 0) && ~isempty( varargin{1} )
                obj.N = length(varargin);
                for ii = 1:obj.N
                    obj.Ts(ii) = size( varargin{ii}, 2 );
                    obj.D = size( varargin{1}, 1 );
                end
                obj.Xdata = zeros( obj.D, sum(obj.Ts) );
                obj.aggTs = [0 cumsum(obj.Ts)];
                for ii = 1:obj.N
                    obj.Xdata( :, obj.aggTs(ii)+1:obj.aggTs(ii+1) ) = varargin{ii};
                end                
            else
                obj.N = 0;
                obj.D = 0;
                obj.aggTs = 0;
                obj.seqNames = {};
            end
            
        end
        
        function retStr = toString(obj )
           retStr = sprintf( 'Emissions Type: %s. Dim: %d', obj.getObsType(), obj.D ); 
        end
        
        function obsTypeStr = getObsType( obj )
            uVals = unique( obj.Xdata );
            if length( uVals ) == 2 && uVals(1) == 0 && uVals(2) == 1
                obsTypeStr = 'Bernoulli';
            elseif all( uVals == int32(uVals) )
                obsTypeStr = 'Multinomial';
            elseif strcmp( class(obj), 'ARSeqData' )
                obsTypeStr = 'AR';
            else
                obsTypeStr = 'Gaussian';
            end
        end
        
        function zTrue = zTrue(obj, ii )
            zTrue = obj.zTrueAll(  obj.aggTs(ii)+1:obj.aggTs(ii+1) );
        end
        
        function Xseq = seq( obj, ii )
            Xseq= obj.Xdata( :, obj.aggTs(ii)+1:obj.aggTs(ii+1) );
        end
        
        function name = name( obj, ii );
            name = obj.seqNames{ii};
        end
        
        function obj = addSeq( obj, seqData, name, zTrue )
            T=size( seqData, 2 );
            obj.N = obj.N+1;
            obj.D = size( seqData,1);
            obj.Ts(end+1) = T;
            obj.aggTs(end+1) = obj.aggTs(end)+T;
            obj.Xdata( :, obj.aggTs(end-1)+1:obj.aggTs(end) ) = seqData;
            if exist( 'name', 'var')
                obj.seqNames{end+1} = name;
            else
                obj.seqNames{end+1} = num2str( obj.N );
            end
            if exist( 'zTrue','var')
                obj.zTrueAll( obj.aggTs(end-1)+1:obj.aggTs(end) )    = zTrue;
            end
        end
    end
end
