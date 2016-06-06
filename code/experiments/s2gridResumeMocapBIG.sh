# WARNING: This script will only use the "ang" group of machines"

Niter=15000

for jobID in 2143090 2143091
do
    for taskID in 1 2 3 4 5
    do

    COMMAND="resumeMocap.sh $jobID $taskID $Niter"

    echo $COMMAND

    #continue;

    qsub $COMMAND

    done

done

echo " "
exit

#######################

runMocapBIG.sh SM+zDD+AnnealLin seq 10000
Your job-array 2143090.1-20:1 ("runMocapBIG.sh") has been submitted

# Restarted as jobs  for first 5 tasks
#  2199353 - 2199357

runMocapBIG.sh SM+zDD+AnnealLin one 10000
Your job-array 2143091.1-20:1 ("runMocapBIG.sh") has been submitted

# Restarted as jobs  for first 5 tasks
#  2199358 - 2199362
