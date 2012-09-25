#include "Eigen/Dense"
#include "mex.h"
#include <math.h>
#include <stdio.h>

typedef Eigen::ArrayXXd MatrixType;
typedef Eigen::ArrayXi StateVectorType;
typedef Eigen::ArrayXd VectorType;
typedef Eigen::ArrayXXi StateMatrixType;


void 
SmoothBack(const MatrixType& transmat, const MatrixType& softev, 
                MatrixType& beta, MatrixType& PM)
{
  int K = (int)softev.rows();
  int T = (int)softev.cols();
  if (beta.cols() != T && beta.rows() != K) 
    beta.resize(K, T);
  if (PM.cols() != T && PM.rows() != K)
    PM.resize(K,T);  
  for (int k = 0; k < K; ++k) 
    beta(k, T-1) = 1.0;
  for (int t = T-2; t >= 0; --t) {
    PM.col(t+1) = beta.col(t+1) * softev.col(t+1);
    beta.col(t) = transmat.matrix() * ( PM.col(t+1) ).matrix();
    // Normalize
    beta.col(t) /= beta.col(t).sum();
  }
  // Get partial marg at first timestep
  PM.col(0) = beta.col(0) * softev.col(0);
  // Debug statements, for reference
  //mexPrintf( "%.3f %.3f beta\n", beta(0,t+1), beta(1,t+1) );
  //mexPrintf( "%.3f %.3f softev\n", softev(0,t+1), softev(1,t+1) );
  //mexPrintf( "%.3f %.3f\n\n", PM(0,t+1), PM(1,t+1) );
}


//
// Returns alpha and loglik
// [alpha ] = function SmoothBackC(transmat, softev)
//
// beta is [K x T]
//
// softev is [K x T]
// transmat is [K x K]
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if (nrhs != 2) {
		mexErrMsgTxt("Needs 2 arguments -- transmat, softev");
		return;
	}

	double* A = mxGetPr(prhs[0]);	
	double* D = mxGetPr(prhs[1]);

	const mwSize* A_dims = mxGetDimensions(prhs[0]);
	const mwSize* D_dims = mxGetDimensions(prhs[1]);

	int K = D_dims[0];	
	int T = D_dims[1];

	if (K != A_dims[0]) {
		mexErrMsgTxt("Softev must be K x T");
		return;
	}


	Eigen::Map<MatrixType> softev(D, K, T);
	Eigen::Map<MatrixType> transmat(A, K, K);
	MatrixType beta = MatrixType::Zero(K,T);
	MatrixType PM   = MatrixType::Zero(K,T);

	SmoothBack(transmat, softev, beta, PM);
	
	double* outputToolPtr;
	plhs[0] = mxCreateDoubleMatrix(K, T, mxREAL);
	outputToolPtr = mxGetPr(plhs[0]);
	memcpy(outputToolPtr, beta.data(), K*T*sizeof(double));

	plhs[1] = mxCreateDoubleMatrix(K, T, mxREAL);
	outputToolPtr = mxGetPr(plhs[1]);
	memcpy(outputToolPtr, PM.data(), K*T*sizeof(double) );

}

