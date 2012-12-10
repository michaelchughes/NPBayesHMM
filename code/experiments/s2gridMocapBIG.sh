# WARNING: This script will only use the "ang" group of machines"

TimeLimit=10000

for initName in 'seq' 'one' #'seq'
do

    for infName in  'SM+zDD+AnnealLin' #'SM+DD'  #'SMnoqrev+DD'
    do

    COMMAND="runMocapBIG.sh $infName $initName $TimeLimit"

    echo $COMMAND

    #continue;

    if [ "$initName" == "cheat" ]
    then
      qsub -t 1-4 $COMMAND
    elif [[ ("$initName" == "one") && ($infName == "Prior") ]]
    then
      continue
    else
      qsub -t 1-20 $COMMAND
    fi

    done

done

echo " "
exit

runMocapBIG.sh SM+zDD+AnnealLin seq 10000
Your job-array 2143090.1-20:1 ("runMocapBIG.sh") has been submitted
runMocapBIG.sh SM+zDD+AnnealLin one 10000
Your job-array 2143091.1-20:1 ("runMocapBIG.sh") has been submitted


####################################
2130597 0.50146 runMocapBI mhughes      qw    11/15/2012 15:26:19                                    1 1-4:1
2130598 0.50135 runMocapBI mhughes      qw    11/15/2012 15:26:19                                    1 1-4:1
2130599 0.50134 runMocapBI mhughes      qw    11/15/2012 15:26:19                                    1 1-4:1
2130600 0.50134 runMocapBI mhughes      qw    11/15/2012 15:26:19                                    1 1-4:1
2130601 0.50133 runMocapBI mhughes      qw    11/15/2012 15:26:19                                    1 1-4:1
2130602 0.50133 runMocapBI mhughes      qw    11/15/2012 15:26:20                                    1 1-4:1
2131131 0.50001 runMocapBI mhughes      qw    11/16/2012 13:41:50                                    1 1-4:1 <<<<<<<<<<  does mem=2GB matter?  NO it doesnot

