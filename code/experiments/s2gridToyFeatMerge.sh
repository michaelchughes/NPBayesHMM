# WARNING: This script will only use the "ang" group of machines"

TimeLimit=3500

for DATANAME in 'AR' #'Gaussian'
do 

for infName in 'Prior' 'SM' 'SMmergeseq'
do

  COMMAND="runToyDataFeatMerge.sh $DATANAME $infName $TimeLimit"

  echo $COMMAND

  #continue;

  qsub -t 1-3 $COMMAND

done
done
exit

##########################################################################  Gauss N=400, T=200
runToyDataFeatMerge.sh Gaussian Prior 3500
Your job-array 1633654.1-5:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh Gaussian SM 3500
Your job-array 1633655.1-5:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh Gaussian SMmergeseq 3500
Your job-array 1633656.1-5:1 ("runToyDataFeatMerge.sh") has been submitted

##########################################################################  AR N=400, T=200
runToyDataFeatMerge.sh AR Prior 3500
Your job-array 1635154.1-3:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh AR SM 3500
Your job-array 1635155.1-3:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh AR SMmergeseq 3500
Your job-array 1635156.1-3:1 ("runToyDataFeatMerge.sh") has been submitted

