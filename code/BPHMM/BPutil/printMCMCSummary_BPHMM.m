function [] = printMCMCSummary_BPHMM( iter, Psi, logPr, algParams )

fprintf( '\t % 5d/%d after %6.0f sec | logPr % .2e | nFeats %4d \n', ...
     iter, algParams.Niter, toc, logPr.all, size(Psi.F,2) );

end