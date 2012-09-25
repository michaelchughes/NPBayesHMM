function logPr = calcMargLogPrData_ARMatrixNormalInvWishart( Xstats, PP )

[D DR K] = size( Xstats.XY );
LOG_PI =  1.144729885849400;

logPr = -Inf( 1, K );
for jj = 1:K
    XX = Xstats.XX(:,:,jj);
    XY = Xstats.XY(:,:,jj);
    YY = Xstats.YY(:,:,jj);
    CC = PP.invAScaleMatrix;
    
    degFreeN = PP.degFree + Xstats.nObs(jj);

    ScaleMatrixN = PP.ScaleMatrix + XX - ( XY/(YY + CC) )*XY';
    %ScaleMatrixN = PP.ScaleMatrix + XX - XY*( (YY+CC)\(XY') );
        
    logPr(jj) = logMvGamma( 0.5*degFreeN, D ) - logMvGamma( 0.5*PP.degFree, D ) ...
                  + 0.5*PP.degFree*log( det( PP.ScaleMatrix ) )   ...
                  - 0.5*degFreeN  *log( det(   ScaleMatrixN ) ) ...
                  + 0.5*D*log( det( CC ) ) ...
                  - 0.5*D*log( det( CC + YY ) );
end

N = sum(Xstats.nObs);
logPr = sum( logPr ) + -0.5*N*D*LOG_PI;
