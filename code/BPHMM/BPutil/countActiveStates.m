function [nActive, F_used, C_used] = countActiveStates( F, stateSeq )

C_used = zeros(size(F));
F_used = zeros(size(F));
for ii = 1:length( stateSeq )
    C_used(ii, :) = histc( stateSeq(ii).z, 1:size(F,2)  );
    F_used(ii, stateSeq(ii).z ) = 1;
end

nActive = sum( sum(F_used,1)   > 0 );