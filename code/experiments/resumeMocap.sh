#!/bin/bash 
#  SYNTAX:  resumeMocap.sh <jobID> <taskID> <Niter>
#$ -S /bin/sh 
#$ -cwd 
# ------ attach job number
#$ -j n
# ------ send to particular queue
#$ -l vlong
# ------ send to particular machine
#$ -q '*@@ang'
# put stdout and stderr files in the right place for your system.
#   NOTE that $TASK_ID is the correct var here
#          but not in rest of script (where SGE_TASK_ID is correct)
#$ -o ../../logs/$JOB_ID.$TASK_ID.out
#$ -e ../../logs/$JOB_ID.$TASK_ID.err

# Edit this line!
PROJHOME=/home/mhughes/git/NPBayesHMM/code/

# ======================================================= SANITIZE INPUT
EXPECTED_ARGS=3
E_BADARGS=65
if [ $# -lt $EXPECTED_ARGS ]  # use -lt instead of <
then
  echo "Usage: `basename $0` <jobID> <taskID> <TimeLimit (seconds)>"
  exit $E_BADARGS
fi
jobID=$1
taskID=$2
Niter=$3


# ======================================================= RUN GRID JOB
matlabpath=/local/projects/matlab/current/bin/
export LD_LIBRARY_PATH=$matlabpath;
#MYPREFIX="dbstop if error; cd $PROJHOME; addpath( genpath( '.'));"
MYPREFIX="cd $PROJHOME; addpath( genpath( '.'));"
MYCMD="ResumeExperiment $jobID $taskID $Niter"

echo "EXECUTING CMD:"
echo "              "$MYPREFIX
echo "              "$MYCMD
$matlabpath/matlab -nojvm -nodesktop -nodisplay -nosplash -r "$MYPREFIX; $MYCMD; exit;"
