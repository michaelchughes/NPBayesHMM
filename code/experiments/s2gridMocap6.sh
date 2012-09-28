# WARNING: This script will only use the "ang" group of machines"

TimeLimit=36000

for initName in 'one' 'unique5'
do

    for infName in 'Prior' 'SM+DD' 'SMnoqrev+DD'
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

