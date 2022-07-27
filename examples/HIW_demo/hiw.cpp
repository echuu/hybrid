#include <RcppArmadillo.h>
#include <cmath>


typedef unsigned int u_int;


// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u);
inline float logmultigamma (unsigned int p, float a);
float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params); 
float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params); 
arma::mat vec2chol(Rcpp::NumericVector& u, Rcpp::List& params);


// [[Rcpp::export]]
inline float psi(Rcpp::NumericVector& u, Rcpp::List& params) {
	return -cov_logprior(u, params) - cov_loglik(u, params); 
}

// [[Rcpp::export]]
float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params) {
	int N = params["N"];
	int D = params["D"];
	arma::mat S = params["S"];
	
	arma::mat Lt = vec2chol(u, params);
	float logdet = 0;
	for (int i = 0; i < D; i++) {
		logdet += std::log(Lt(i,i));
	}
	float loglik = -0.5 * N * D * std::log(2*M_PI) + N * logdet - 
					0.5 * arma::trace(Lt.t() * Lt * S);
	return loglik;
}

// [[Rcpp::export]]
arma::mat vec2chol(Rcpp::NumericVector& u, Rcpp::List& params) {

	u_int D      = params["D"];
	arma::mat G  = params["G"];

	arma::mat Lt(D, D, arma::fill::zeros);
	u_int p = 0;
	for (u_int c = 0; c < D; c++) {
		for (u_int r = 0; r <= c; r++) {
			if (G(r,c) > 0) {
				Lt(r, c) = u[p++];	
			}
		}
	}
	return Lt;
}



// [[Rcpp::export]]
Rcpp::NumericVector chol2vec(arma::mat& M, Rcpp::List& params) {

	u_int D      = params["D"];
	u_int D_u    = params["D_u"];
	arma::mat G  = params["G"];
	u_int k      = 0; // used to increment the index in the vector u 
	Rcpp::NumericVector u(D_u);
	for (u_int c = 0; c < D; c++) {
		for (u_int r = 0; r <= c; r++) {
			if (G(r,c) > 0) {
				u[k++] = M(r, c);
			}
		}
	}
	return u;
}


// [[Rcpp::export]]
Rcpp::NumericVector grad(Rcpp::NumericVector& u, Rcpp::List& params) {
	arma::vec nu     = params["nu"];
	arma::vec xi     = params["xi"];
	arma::mat V      = params["V"];
	arma::mat S      = params["S"]; 
	u_int D          = params["D"];   // dimension of the cholesky factor
	u_int D_u        = params["D_u"]; // dimension of the parameter
	u_int N          = params["N"];

	arma::mat Lt = vec2chol(u, params);

	arma::mat diag_terms(D, D, arma::fill::zeros);
	for (u_int i = 0; i < D; i++) { 
		diag_terms(i,i) = (xi[i] + N) / Lt(i,i); 
	}
	
	arma::mat grad_mat(D, D, arma::fill::zeros);
	grad_mat = - diag_terms + Lt * S + Lt;

	Rcpp::NumericVector grad_vec = chol2vec(grad_mat, params);
	return grad_vec;
	// return Lt;
}


// [[Rcpp::export]]
arma::mat hess(Rcpp::NumericVector& u, Rcpp::List& params) {

	arma::vec nu      = params["nu"];
	arma::vec xi      = params["xi"];
	arma::mat V       = params["V"];
	arma::mat S       = params["S"]; 
	unsigned int D    = params["D"];
	unsigned int D_u  = params["D_u"];
	u_int N           = params["N"];
	arma::mat ind_mat = params["t_ind"];

	// initialize the output matrix
	arma::mat H(D_u, D_u, arma::fill::zeros);

	// reconstruct upper cholesky factor -- to be moved into another function
	arma::mat Lt = vec2chol(u, params);


	unsigned int i, j, k, l, r, c;
	// arma::mat test(D_u, 2, arma::fill::zeros);

	for (r = 0; r < D_u; r++) {
		i = ind_mat(r, 0) - 1; // subtract one to account for 0-index
		j = ind_mat(r, 1) - 1;
		c = r;
		while (c < D_u) {

			k = ind_mat(c, 0) - 1; // row of 2nd order partial
			l = ind_mat(c, 1) - 1; // col of 2nd order partial

			if (i != k) {
				H(r,c) = 0;
				H(c,r) = 0;
			} else if (i == j && k == i && l > j) {
				H(r,c) = -S(l,j);
				H(c,r) = H(r,c);
			} else if (i == j && j == k && k == l) {
				// -1/Lt[i,i]^2 * (N + xi[i]) - S[i,i] - 1
				H(r,c) = -(xi(i) + N) / std::pow(Lt(i,i), 2) - S(i,i) - 1;
				H(c,r) = H(r,c);
			} else if (i != j && k == i && l == j) {
				H(r,c) = -S(l, j) - 1;
				H(c,r) = H(r,c);
			} else if (i != j && k == i && l > j) {
				H(r,c) = -S(l, j);
				H(c,r) = H(r,c);
			}
			c++;
		} // end inner while loop
		// return ind_mat;

	} // end for() populating hessian matrix

	// H = 0.5 * (H + H.t());

	return -H;
}


// [[Rcpp::export]]
float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params) {

	arma::mat V     = params["V"];
	unsigned int b  = params["b"];
	unsigned int D  = params["D"];
	arma::vec nu    = params["nu"];
	arma::mat G  = params["G"];
	// arma::mat Lt(D, D, arma::fill::zeros);

	arma::mat Lt = vec2chol(u, params);
	float logprior = 0;
	// TODO: compute upper diagonal part of the log prior
	for (int r = 0; r < D; r++) {
		for (int c = (r+1); c < D; c++) {
			if (G(r,c) > 0) {
				logprior += -0.5 * std::log(2 * M_PI) - 0.5 * std::pow(Lt(r,c), 2);
			}
		}
	}


	// compute diagonal part of the log prior
	for (unsigned int i = 0; i < D; i++) {
		logprior += -0.5*(b + nu(i)) * std::log(2) - 
						lgamma(0.5*(b + nu(i))) +  
						(b + nu(i) - 2) * std::log(Lt(i,i)) - 
						0.5 * std::pow(Lt(i,i), 2) + std::log(2 * Lt(i,i));
	}

	return logprior;
}

/* Populars lower diagonal of N x N matrix L with values from u,
	params:
		L: N x N arma matrix
		u: numeric vector
	returns:
		sum of log diag values
*/
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u) {
	if (L.n_rows != L.n_cols)
		throw std::invalid_argument("Unequal dimensions");

	unsigned int D = L.n_rows;
	unsigned int k = 0;	
	float logDiagL = 1;

	for (unsigned int j = 0; j < D; j++) {
		logDiagL *= u[k];
		for (unsigned i = j; i < D; i++) {
			L(i, j) = u[k++];	
		}
	}
	return std::log(logDiagL);
}

/* log multivariate gamma function
	params:
*/
inline float logmultigamma (unsigned int p, float a) {
	float f = 0.25 * p * (p - 1) * std::log(M_PI);
	for (unsigned int i = 1; i <= p; i++)
		f += lgamma(a+0.5-(0.5 * i));
	return f;
}

