// informME: An information-theoretic pipeline for WGBS data
// Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// or see <http://www.gnu.org/licenses/>.
//
// computeMCtransProbs.cpp

#include <sstream>
#include <iostream>
#include <Eigen/MPRealSupport>
#include "mpreal.h"
#include "mex.h"

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Computation of the initial and transition probabilities 
	 * of the alternative inhomogeneous Markov chain  
	 * representation of the 1D Ising model
     	 *
     	 * Author: W. Garrett Jenkinson
	 * 
     	 * Last Modified: November 18, 2016
     	 *
     	 * usage:
	 * [p1,transProbs] = computeMCtransProbs(An,Cn,logZ1,logZ0,logZ);
     	 *
     	 * INPUTS:
     	 * An is a Nx1 vector of A_n parameters for n=1,...,N
	 * Cn is a (N-1)x1 vector of C_n parameters for n=2,...,N
	 * z1 is a Nx1 vector of values Z_n(x_n) for x_n=1 and n=1,...,N
     	 * z0 is a Nx1 vector of values Z_n(x_n) for x_n=0 and n=1,...,N
     	 * Z is the partition function
    	 * 
     	 * OUTPUTS:
     	 * p1 is the marginal probability P(X_1=0)
     	 * transProbs is a (N-1)x2 matrix
	 * - First  column has probability P(X_{n+1}=0|X_n=0) for n=1,2,...,N-1
	 * - Second column has probability P(X_{n+1}=0|X_n=1) for n=1,2,...,N-1
	 *
	 * Note that the Ising distribution can be calculated from these probabilities: 
	 * Pr(X_1,...X_N) = Pr(X_1)\prod_{n=1}^{N-1} Pr(X_{n+1}|X_n)
	 * 
	 * More importantly these probabilities can be used to iteratively draw an
	 * exact sample from the Ising distribution without using Markov Chain 
	 * Monte Carlo (MCMC) sampling.
	 */
	 
    // Set Multiple Precision Code
    using mpfr::mpreal;
    
    // Required precision of computations in decimal digits
    // Play with it to check different precisions
    const int digits = 200;
    
    // Initialize Multiple Precision Number 
    const mpreal zero = 0.0;
    
    // Set default precision for all subsequent computations
    // MPFR accepts precision in bits - so do the conversion
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    
    // Declare matrix and vector types with multi-precision scalar type
    typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
    typedef Matrix<mpreal,Dynamic,1>        VectorXmp;
    
    // Initialize matlab input/output vars
	
    int N = mxGetNumberOfElements(prhs[0]); //number of elements of first rhs argument
    const double *An  = (double *)mxGetData(prhs[0]);  // read in first rhs argument
    const double *Cn = (double *)mxGetData(prhs[1]);  // read in second rhs argument
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); // create first output argument in matlab
    plhs[1] = mxCreateDoubleMatrix(N-1,2,mxREAL); // create second output argument in matlab
    
    // Make mpreal outputs (later converted to strings when passed to matlab)
    MatrixXmp transProbs(N-1,2);
    mpreal p1;

    // Make logZ0 and logZ1 vectors 
    VectorXmp logZ0(N);
    VectorXmp logZ1(N);
    
    // Make logZ scalar
    mpreal logZ;
        
    // Begin computation of logZ's
 
    // Calculate first boundary values
    logZ1(N-1)=0.0;
    logZ0(N-1)=0.0;
    
    // Recurse through non-boundary values
    for (int n=N-1;n>1;n--){
        logZ1(n-1) = log(exp(-An[n]-Cn[n-1]+logZ0(n))
                         + exp(An[n]+Cn[n-1]+logZ1(n)));
        logZ0(n-1) = log(exp(-An[n]+Cn[n-1]+logZ0(n))
                        + exp(An[n]-Cn[n-1]+logZ1(n)));
    }
    
    // Calculate last boundary values
    logZ1(0) = log(exp(An[0]-An[1]-Cn[0]+logZ0(1))
                  + exp(An[0]+An[1]+Cn[0]+logZ1(1)));
    logZ0(0) = log(exp(-An[0]-An[1]+Cn[0]+logZ0(1))
                  + exp(-An[0]+An[1]-Cn[0]+logZ1(1)));
    
    // Compute log(Z_1(0)+Z_1(1))
    logZ = log(exp(logZ0(0)) + exp(logZ1(0)));
    
    // Compute boundary condition 
    p1 = exp(logZ0(0)-logZ); // P(X_N=0)
    transProbs(0,0) = exp(-An[0]-An[1]+Cn[0]+logZ0(1)-logZ0(0));     
    transProbs(0,1) = exp(An[0]-An[1]-Cn[0]+logZ0(1)-logZ1(0));
     
    // Use backwards recursion to compute transition probabilities
    for (int n=2;n<N;n++){
        transProbs(n-1,0) = exp(-An[n]+Cn[n-1]+logZ0(n)-logZ0(n-1));         
        transProbs(n-1,1) = exp(-An[n]-Cn[n-1]+logZ0(n)-logZ1(n-1));
    }

    // Write output to Matlab
    double *outputProb;
    outputProb = mxGetPr(plhs[0]);
    outputProb[0]=p1.toDouble();
    
    double *outputMatrix;
    outputMatrix = mxGetPr(plhs[1]);
    for (int n = 0; n < N-1; n++) {
        outputMatrix[n]     = transProbs(n,0).toDouble();
        outputMatrix[n+N-1] = transProbs(n,1).toDouble();

    }

}
