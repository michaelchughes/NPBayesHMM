# WARNING: This script will only use the "ang" group of machines"

TimeLimit=3500
DATANAME='Gaussian'

for infName in 'Prior' 'DD' 'SM'
do

  COMMAND="runToyDataFeatCreation.sh $DATANAME $infName $TimeLimit"

  echo $COMMAND

  #continue;

  qsub -t 1-5 $COMMAND

done

exit


#
runToyDataFeatCreation.sh Gaussian Prior 3500
Your job-array 1511213.1-5:1 ("runToyDataFeatCreation.sh") has been submitted
runToyDataFeatCreation.sh Gaussian DD 3500
Your job-array 1511214.1-5:1 ("runToyDataFeatCreation.sh") has been submitted
runToyDataFeatCreation.sh Gaussian SM 3500
Your job-array 1511215.1-5:1 ("runToyDataFeatCreation.sh") has been submitted

