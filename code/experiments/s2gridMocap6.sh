# WARNING: This script will only use the "ang" group of machines"

TimeLimit=36000

for initName in 'one' #'unique5'
do

    for infName in 'SM+zDD' 'SM+DD+Anneal' 'SM+zDD+Anneal' #'Prior' 'SM+DD' 'SMnoqrev+DD'
    do

    COMMAND="runMocap6.sh $infName $initName $TimeLimit"

    echo $COMMAND

    #continue;

    if [ "$initName" == "cheat" ]
    then
      qsub -t 1-10 $COMMAND
    elif [[ ("$initName" == "one") && ($infName == "Prior") ]]
    then
      continue
      #qsub -t 1 $COMMAND    
    else
      qsub -t 1-10 $COMMAND
    fi

    done

done

echo " "
exit

runMocap6.sh SM+zDD one 36000
Your job-array 1876678.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh SM+DD+Anneal one 36000
Your job-array 1876679.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh SM+zDD+Anneal one 36000
Your job-array 1876680.1-10:1 ("runMocap6.sh") has been submitted


Skipped Prior one since it is worthless.
runMocap6.sh SM+DD one 36000
Your job-array 1515157.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh SMnoqrev+DD one 36000
Your job-array 1515158.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh Prior unique5 36000
Your job-array 1515159.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh SM+DD unique5 36000
Your job-array 1515160.1-10:1 ("runMocap6.sh") has been submitted
runMocap6.sh SMnoqrev+DD unique5 36000
Your job-array 1515161.1-10:1 ("runMocap6.sh") has been submitted

