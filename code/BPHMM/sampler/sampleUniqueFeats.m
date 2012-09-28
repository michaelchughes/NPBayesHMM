function [Psi, Stats, RhoTerms] = sampleUniqueFeats( Psi, data, algParams )
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

Stats.ADD.nAccept = 0;
Stats.ADD.nTotal = 0;
Stats.DEL.nAccept = 0;
Stats.DEL.nTotal = 0;

% Fill in new values!!
for ii = 1:size( Psi.F, 1)
    [Psi, RhoTerms] = sampleSingleFeatEntry_UniqueRJ( ii, Psi, data, algParams );
    
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

% ===================================================== Record Descr Stats
% If MoveType=0 (death), then if accepted fNew should be 0
%  likewise, MoveType=1(birth) should if accepted make fNew=1
% MoveStats.nAccept = sum( fNew == MoveType01 );
% MoveStats.nTrial  = nObj;
% MoveStats.BIRTH.nAccept = sum( fNew( MoveType01 == 1 ) == 1 );
% MoveStats.BIRTH.nTrial  = sum( MoveType01 == 1 );
% MoveStats.DEATH.nAccept = sum( fNew( MoveType01 == 0 ) == 0 );
% MoveStats.DEATH.nTrial  = sum( MoveType01 == 0 );
