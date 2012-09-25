if ismac
   mex -v -I/usr/local/include/eigen3/ SampleHMMStateSeqWithQsC.cpp
   mex -v -I/usr/local/include/eigen3/ SampleHMMStateSeqC.cpp   
   mex -v -I/usr/local/include/eigen3/ MySmoothBackC.cpp   
   mex -v -I/usr/local/include/eigen3/ FilterFwdC.cpp   
else
   mex -v -I/home/mhughes/code/Eigen/include/eigen3/ SampleHMMStateSeqWithQsC.cpp
   mex -v -I/home/mhughes/code/Eigen/include/eigen3/ SampleHMMStateSeqC.cpp   
   mex -v -I/home/mhughes/code/Eigen/include/eigen3/ MySmoothBackC.cpp
   mex -v -I/home/mhughes/code/Eigen/include/eigen3/ FilterFwdC.cpp
end
