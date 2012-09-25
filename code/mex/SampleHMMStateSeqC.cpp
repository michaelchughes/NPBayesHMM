#include "Eigen/Dense"
#include "mex.h"
#include <math.h>
#include "mersenneTwister2002.c"

typedef Eigen::ArrayXXd MatrixType;
typedef Eigen::ArrayXi StateVectorType;
typedef Eigen::ArrayXd VectorType;
typedef Eigen::ArrayXXi StateMatrixType;

int draw_discrete( VectorType& ps )
{
    double total = ps.sum();
    int K = (int)ps.size();
    
    double r = total*genrand_double();
    double cursum = ps(0);
    int newk = 0;
    while ( r >= cursum && newk < K-1) {
        newk++;
        cursum += ps[newk];
    }
    if ( newk < 0 || newk >= K ) {
        mexErrMsgTxt( "   somethings not right here.");
        return -1;
    }
    return newk;
}

void 
SampleZ( const MatrixType& transmat, const Eigen::Map<VectorType> pi_init, 
        const MatrixType& PM, MatrixType& Z)
{
    int K = (int)PM.rows();
    int T = (int)PM.cols();
    VectorType ps = VectorType::Zero(K);
    ps = pi_init;
    ps *= PM.col(0);
    Z(0) = draw_discrete( ps );
    for (int t = 1; t<T; t++) {
        ps = PM.col(t);
        ps *= transmat.row(Z(t-1));
        Z(t) = draw_discrete( ps );
    }
    Z += 1;
}


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
	if (nrhs != 4) {
		mexErrMsgTxt("Needs 4 arguments -- transmat, softev, init, SEED");
		return;
	}

	double* A = mxGetPr(prhs[0]);	
	double* D = mxGetPr(prhs[1]);
    double* Q = mxGetPr(prhs[2]);
    int SEED = (int)*mxGetPr(prhs[3]);
    
    // Initialize rand generator with the SEED ---------------------------
    // int SEED = (int) mxGetScalar( SEED_IN );
    init_genrand( SEED );
    
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
    Eigen::Map<VectorType> init(Q, K );
	MatrixType beta = MatrixType::Zero(K,T);
	MatrixType PM   = MatrixType::Zero(K,T);
    MatrixType Z    = MatrixType::Zero(1,T);
    
	SmoothBack(transmat, softev, beta, PM);
	SampleZ( transmat, init, PM, Z);
    
    
	double* outputToolPtr;
	plhs[0] = mxCreateDoubleMatrix(1, T, mxREAL);
	outputToolPtr = mxGetPr(plhs[0]);
	memcpy(outputToolPtr, Z.data(), T*sizeof(double));

}

