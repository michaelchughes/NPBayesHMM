function [Psi, Stats, RhoTerms] = sampleUniqueFeats( Psi, data, algParams, doZMove, objIDs )
% Sample *unique* features for each time series, 
%   using reversible jump moves that propose adding/deleting unique feat.
%   from each sequence.
% The proposal distribution for reversible jump is defined in "algParams'
%   can be any of:
%      'prior'  : emit param thetaStar drawn from prior
%      'data-driven' : emit param thetaStar draw from data posterior
%OUTPUT
%  Psi : resulting model config
%  Stats : struct that summarizes reversible jump performance
%           counts # birth (ADD) and death (DEL) attempts and acceptances

if ~exist( 'objIDs', 'var')
    objIDs = 1:data.N;
end


Stats.ADD.nAccept = 0;
Stats.ADD.nTotal = 0;
Stats.DEL.nAccept = 0;
Stats.DEL.nTotal = 0;

for ii = objIDs
    if doZMove
        [Psi, RhoTerms] = sampleSingleFeat_UniqueRJStateSeq( ii, Psi, data, algParams );
    else
        [Psi, RhoTerms] = sampleSingleFeatEntry_UniqueRJ( ii, Psi, data, algParams );
    end
    
    if RhoTerms.doBirth
        Stats.ADD.nTotal = Stats.ADD.nTotal+1;
        if RhoTerms.doAccept;
            Stats.ADD.nAccept = Stats.ADD.nAccept+1;
        end
        
    else
        Stats.DEL.nTotal = Stats.DEL.nTotal+1;
        if RhoTerms.doAccept;
            Stats.DEL.nAccept = Stats.DEL.nAccept+1;
        end
    end
end
