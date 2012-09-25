#!/bin/bash 
# CompileMEX.sh
# Run this script ONCE (when you first clone the git repository)
#  to compile the fast MEX routines for HMM dynamic programming.

if [ -e "code/EigenLibrary.path" ]
then
  read eigenLibPath < code/EigenLibrary.path
else
  read -p "Enter the path to Eigen library:" eigenLibPath  
  echo $eigenLibPath > code/EigenLibrary.path
fi

if [ ! -d "$eigenLibPath" ]
then
  echo "ERROR: Bad path in code/EigenLibrary.path."
  echo "HINT: good path is a directory of the form  </path/to/Eigen>/include/eigen3/"
  echo "   where </path/to/Eigen> will vary by user"
  echo "Please fix the file and run this script again to compile."
  exit
fi

# =======================================================================  MAKE MEX BINARIES
echo "Compiling MEX binaries... "
echo " using Eigen path:" $eigenLibPath
cd code/mex/
mex -I$eigenLibPath SampleHMMStateSeqWithQsC.cpp
mex -I$eigenLibPath SampleHMMStateSeqC.cpp
mex -I$eigenLibPath MySmoothBackC.cpp
mex -I$eigenLibPath FilterFwdC.cpp
cd ../../

echo "... successfully compiled fast MEX routines for HMM dynamic programming."
