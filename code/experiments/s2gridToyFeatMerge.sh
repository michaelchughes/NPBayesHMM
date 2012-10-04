#!/bin/sh
# s2gridToyFeatMerge.sh <dataName> <doSubmit>
# WARNING: This script will only use the "ang" group of machines"

EXPECTED_ARGS=1
E_BADARGS=65
if [ $# -lt $EXPECTED_ARGS ]  # use -lt instead of <
then
  echo "Usage: `basename $0` <dataName=AR>"
  exit $E_BADARGS
fi

TimeLimit=3500
DATANAME=$1

for infName in 'Prior' 'SM' 'SMmergeseq'
do

  COMMAND="runToyDataFeatMerge.sh $DATANAME $infName $TimeLimit"

  echo $COMMAND

  if [ $# -eq 2 ]
  then
    qsub -t 1-$2 $COMMAND
  else
    continue
  fi

done
exit
##########################################################################  AR N=200, T=1000, kappa=50 (nearly same as orig. experiment)
runToyDataFeatMerge.sh AR Prior 3500
Your job-array 1637996.1-2:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh AR SM 3500
Your job-array 1637997.1-2:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh AR SMmergeseq 3500
Your job-array 1637998.1-2:1 ("runToyDataFeatMerge.sh") has been submitted


##########################################################################  Gauss N=600, T=200
runToyDataFeatMerge.sh Gaussian Prior 3500
Your job-array 1637993.1-2:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh Gaussian SM 3500
Your job-array 1637994.1-2:1 ("runToyDataFeatMerge.sh") has been submitted
runToyDataFeatMerge.sh Gaussian SMmergeseq 3500
Your job-array 1637995.1-2:1 ("runToyDataFeatMerge.sh") has been submitted

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

OLD REPEAT JOB        jobIDs = [806212:806214];
