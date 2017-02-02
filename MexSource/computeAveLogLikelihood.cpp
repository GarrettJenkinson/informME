//   informME: An information-theoretic pipeline for WGBS data
//   Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software Foundation,
//   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
//   or see <http://www.gnu.org/licenses/>.
//
// computeAveLogLikelihood.cpp

#include "mex.h"
#include "mpreal.h"
#include <iostream>
#include <Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Computation of the average "marginalized" log-likelihood of 
	 * a set of methylation observations for a given Ising model
     	 *
     	 * Author: W. Garrett Jenkinson
	 * 
     	 * Last Modified: November 18, 2016
     	 *
     	 * usage:
	 * avelogLikelihood = computeAveLogLikelihood(An,Cn,dataMatrix',CpGstart,CpGend);
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
    int N = mxGetNumberOfElements(prhs[0]); // number of elements of first rhs argument
    int K = mxGetNumberOfElements(prhs[3]); // number of elements of fourth rhs argument
    const double *An = (double *)mxGetData(prhs[0]); // read in first rhs argument
    const double *Cn = (double *)mxGetData(prhs[1]); // read in second rhs argument
    const int    *dataMatrix   = (int *)mxGetData(prhs[2]); // read in third rhs argument
    const int    *CpGstart     = (int *)mxGetData(prhs[3]); // read in fourth rhs argument
    const int    *CpGend       = (int *)mxGetData(prhs[4]); // read in fifth rhs argument
    
    plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL); // create first output argument in matlab
    double *avelogLikelihood = (double *)mxGetData(plhs[0]);    // bind this first output argument in C++
    
    // Make logZ0 and logZ1 vectors 
    VectorXmp logZ1(N);
    VectorXmp logZ0(N);
    VectorXmp logZ1tilde(N);
    VectorXmp logZ0tilde(N);
    
    // Make scalars
    mpreal logZ;
    mpreal logMargLikeTemp=0.0;
    mpreal logAveLikeTemp =0.0;
    int r,s;
    int x_r_rPLUSs[N];
    
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
                  + exp( An[0]+An[1]+Cn[0]+logZ1(1)));
    logZ0(0) = log(exp(-An[0]-An[1]+Cn[0]+logZ0(1))
                  + exp(-An[0]+An[1]-Cn[0]+logZ1(1)));
    
    // Compute log(Z_1(0)+Z_1(1))
    logZ = log(exp(logZ0(0)) + exp(logZ1(0)));
    
    // Begin computation of logZtilde's
        
    // Calculate first two boundary values
    logZ1tilde(0) = 0.0;
    logZ0tilde(0) = 0.0;

    logZ1tilde(1) = log(exp(-An[0]+An[1]-Cn[0]+zero) 
                            +exp(An[0]+An[1]+Cn[0]+zero));
    logZ0tilde(1) = log(exp(-An[0]-An[1]+Cn[0]+zero) 
                            +exp(An[0]-An[1]-Cn[0]+zero));
    
    // Recurse through non-boundary values
    for (int n=3;n<N+1;n++){
        logZ1tilde(n-1) = log(exp( An[n-1]-Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(An[n-1]+Cn[n-2]+logZ1tilde(n-2)));
        logZ0tilde(n-1) = log( exp(-An[n-1]+Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(-An[n-1]-Cn[n-2]+logZ1tilde(n-2)));
    }

    // Begin computing likelihoods for each observation
    for (int k=0;k<K;k++){

        r = CpGstart[k];
        s = CpGend[k]-CpGstart[k];
      
        for (int CpG = r;CpG<r+s+1;CpG++){
             
            x_r_rPLUSs[CpG-r]=dataMatrix[CpG-1+k*N]; // dataMatrix is transposed before input
                                                     // so the matrix is stored linearly by rows
        }        
        
        // Set to log(1/Z)
        logMargLikeTemp = -logZ;

        // Add log[\tilde{Z}_{r}(x_r)]
        if (x_r_rPLUSs[0]>0){ //X(r)=1
            logMargLikeTemp = logMargLikeTemp + logZ1tilde(r-1);
        }else {              //X(r)=0
            logMargLikeTemp = logMargLikeTemp + logZ0tilde(r-1);
        }

        // Add log[Z_{r+s}(x_{r+s})]
        if (x_r_rPLUSs[s]>0){ //X(r+s)=1
            logMargLikeTemp = logMargLikeTemp + logZ1(r+s-1);
        }else {               //X(r+s)=0
            logMargLikeTemp = logMargLikeTemp + logZ0(r+s-1);
        }

        // Add by \sum_{n=r}^{r+s-1} log[\phi(x_n,x_{n+1})]
        if (s>0){ // otherwise this product is empty
            if (r>1){ // no boundary \phi term needed
                for (int iter=1;iter<s+1;iter++){ 
                    logMargLikeTemp=logMargLikeTemp + (An[iter+r-1]+Cn[iter+r-2]*(2*x_r_rPLUSs[iter-1]-1))*(2*x_r_rPLUSs[iter]-1);
                }

            }else{ // must compute boundary \phi term since r=1

                // n=1;
                logMargLikeTemp=logMargLikeTemp + (2*x_r_rPLUSs[0]-1)*An[0]
                                +(2*x_r_rPLUSs[1]-1)*(An[1] +Cn[0]*(2*x_r_rPLUSs[0]-1));

                for (int iter2=2;iter2<s+1;iter2++){ 
                    logMargLikeTemp=logMargLikeTemp + (An[iter2]+Cn[iter2-1]*(2*x_r_rPLUSs[iter2-1]-1) )*(2*x_r_rPLUSs[iter2]-1);

                }

            }

        }
        
        logAveLikeTemp = logAveLikeTemp + logMargLikeTemp;
   
    }
    
    logAveLikeTemp = logAveLikeTemp/K;
    avelogLikelihood[0] = logAveLikeTemp.toDouble(); 

}

