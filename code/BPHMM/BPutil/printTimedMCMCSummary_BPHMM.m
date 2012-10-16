function [] = printTimedMCMCSummary_BPHMM( iter, Psi, logPr, algParams )

nDigits = length( num2str(algParams.TimeLimit) );

FMT_STR = sprintf( '%d', nDigits );
FMT_STR = [ '%' FMT_STR '.0f/%' FMT_STR '.0f' ];

if isfield( Psi, 'invTemp' )
fprintf( ['\t ' FMT_STR ' sec | iter % 5d | logPr % .2e | nFeats %4d | 1/temp %.2f\n'], ...
     toc, algParams.TimeLimit, iter, logPr.all, size(Psi.F,2), Psi.invTemp );
else
fprintf( ['\t ' FMT_STR ' sec | iter % 5d | logPr % .2e | nFeats %4d \n'], ...
     toc, algParams.TimeLimit, iter, logPr.all, size(Psi.F,2) );
end

end