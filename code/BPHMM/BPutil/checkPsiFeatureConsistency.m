function boolresult = checkPsiFeatureConsistency( Psi )
boolresult = false;
for ii = 1:size( Psi.F,1)
       kii = unique( Psi.stateSeq(ii).z );
       for kk = kii
           assert( Psi.F(ii,kk) == 1, 'Bad state sequence assignment!' );
       end
       assert( all( find(Psi.F(ii,:) ) == Psi.TransM.seq(ii).availFeatIDs ), 'Bad feature translation!' );
end
boolresult = true;