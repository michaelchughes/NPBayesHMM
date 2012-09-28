function qs = buildRJMoveDistr( uCur, Kii, PrBirth )
% Produce discrete distribution over move options for reversible jump
%   Note that given current state of uCur unique feats,
%   options are birth, delete feat 1, delete feat 2, ... delete feat uCur
%   If uCur=1, we give birth with prob. one
% OUTPUT
%  qs = 1 x (1 + uCur ) row vector
%   where
%      qs( 1 + k ) = prob. of killing unique feature k
%                       e.g. unique_features_ii(k)
%      qs(1) = prob. of giving birth to new feature

% if we only have one feature, never take it away!
if Kii==1
   qs = zeros(1, uCur);
   qs(1) = 1;
else    
    qs = ( 1 - PrBirth( uCur ) )/uCur * ones(1, uCur+1 );
    qs(1) = PrBirth( uCur );
    qs = qs./sum(qs);
end
