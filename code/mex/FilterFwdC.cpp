#include "Eigen/Dense"
#include "mex.h"
#include <math.h>

typedef Eigen::ArrayXXd MatrixType;
typedef Eigen::ArrayXi StateVectorType;
typedef Eigen::ArrayXd VectorType;
typedef Eigen::ArrayXXi StateMatrixType;


void 
FilterFwd(const Eigen::Map<MatrixType>& transmat, const Eigen::Map<MatrixType>& softev, 
          const Eigen::Map<VectorType>& init, double& loglik, MatrixType& alpha)
{
    int T = (int) softev.cols();
    int K = (int) softev.rows();

    if (alpha.cols() != T && alpha.rows() != K) {
        alpha.resize(K, T);
    }
    VectorType scale = VectorType::Zero(T);
    Eigen::MatrixXd at = transmat.matrix().transpose();

    alpha.col(0) = init * softev.col(0);
    scale(0) = alpha.col(0).sum();
    alpha.col(0) /= scale(0);

    for (int t = 1; t < T; ++t) {
        alpha.col(t) = (at.matrix() * alpha.col(t-1).matrix()).array();
        alpha.col(t) *= softev.col(t);
        scale(t) = alpha.col(t).sum();
        alpha.col(t) /= scale(t);
    }
    loglik = scale.log().sum();
}

//
// Returns alpha and loglik
// [alpha, loglik] = function FilterFwdC(transmat, softev, init)
//
// alpha is [K x T]
// loglik is [1x1]
//
// softev is [K x T]
// transmat is [K x K]
// init is [K x 1]  
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if (nrhs != 3) {
		mexErrMsgTxt("Needs 3 arguments -- transmat, softev, init");
		return;
	}

	double* A = mxGetPr(prhs[0]);	
	double* D = mxGetPr(prhs[1]);
	double* Pi = mxGetPr(prhs[2]);

	const mwSize* A_dims = mxGetDimensions(prhs[0]);
	const mwSize* D_dims = mxGetDimensions(prhs[1]);

	int K = D_dims[0];	
	int T = D_dims[1];

	if (K != A_dims[0]) {
		mexErrMsgTxt("Softev must be K x T");
		return;
	}

	double loglik;

	Eigen::Map<MatrixType> softev(D, K, T);
	Eigen::Map<VectorType> init(Pi, K);
	Eigen::Map<MatrixType> transmat(A, K, K);
	MatrixType alpha = MatrixType::Zero(K,T);

	FilterFwd(transmat, softev, init, loglik, alpha);
	
	double* outputToolPtr;
	plhs[0] = mxCreateDoubleMatrix(K, T, mxREAL);
	outputToolPtr = mxGetPr(plhs[0]);
	memcpy(outputToolPtr, alpha.data(), K*T*sizeof(double));

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	outputToolPtr = mxGetPr(plhs[1]);
	outputToolPtr[0] = loglik;
}
