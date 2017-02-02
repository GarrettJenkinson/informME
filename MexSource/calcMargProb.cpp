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
// calcMargProb.cpp

#include <iostream>
#include <Eigen/MPRealSupport>
#include "mpreal.h"
#include "mex.h"

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Computation of the joint probability of a contiguous set 
	 * of CpG sites
	 *
	 * Author: W. Garrett Jenkinson
	 *
	 * Last Modified: November 18, 2016
	 *
	 * usage:
	 * logMargProb = calcMargProb(r,s,x_r_rPLUSs,logZ1,logZ0,
	 *                       logZ,logZ1tilde,logZ0tilde,An,Cn)
     	 *
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
    
    // Initialize matlab input output vars
    
    int N = mxGetNumberOfElements(prhs[8]); // number of elements of nineth rhs argument
    const int *r               = (int *)mxGetData(prhs[0]); // read in first rhs argument
    const int *s               = (int *)mxGetData(prhs[1]); // read in second rhs argument
    const int *x_r_rPLUSs      = (int *)mxGetData(prhs[2]); // read in third rhs argument
    const double *An  = (double *)mxGetData(prhs[8]); // read in nineth rhs argument
    const double *Cn = (double *)mxGetData(prhs[9]); // read in tenth rhs argument
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); // create first output argument in matlab
    double *logMargProbDoub = mxGetPr(plhs[0]); // bind this first output argument in C++
    
    int n;
    
    // Compute logZs
    
    // Make logZ0 and logZ1 vectors 
    VectorXmp logZ1(N);
    VectorXmp logZ0(N);
    
    // Make logZ scalar
    mpreal logZ;
    mpreal logMargProb;

    // Begin computation of logZ's

    // Calculate first boundary values
    logZ1(N-1)=0.0;
    logZ0(N-1)=0.0;
    
    // Recurse through non-boundary values
    for (int n=N-1;n>1;n--){
        logZ1(n-1) = log(exp(-An[n]-Cn[n-1]+logZ0(n))
                        + exp(An[n]+Cn[n-1]+logZ1(n)));
        
        logZ0(n-1) = log(exp(-An[n]+Cn[n-1]+logZ0(n))
                        + exp( An[n]-Cn[n-1]+logZ1(n)));
    }
    
    // Calculate last boundary values
    
    logZ1(0) = log(exp(An[0]-An[1]-Cn[0]+logZ0(1))
                  + exp(An[0]+An[1]+Cn[0]+logZ1(1)));
    
    logZ0(0) = log(exp(-An[0]-An[1]+Cn[0]+logZ0(1))
                  + exp(-An[0]+An[1]-Cn[0]+logZ1(1)));
    
    // Compute log(Z_1(0)+Z_1(1))
    logZ = log(exp(logZ0(0)) + exp(logZ1(0)));
    
    // Compute logZtildes
      
    // Make logZ0tilde and logZ1tilde vectors 
    VectorXmp logZ0tilde(N);
    VectorXmp logZ1tilde(N);
    
    // Make logZ scalar
    mpreal logZtilde;

    // Calculate first two boundary values
    logZ1tilde(0) = 0.0;
    logZ0tilde(0) = 0.0;
            
    logZ1tilde(1) = log(exp(-An[0]+An[1]-Cn[0]+zero) 
                            +exp(An[0]+An[1]+Cn[0]+zero));
   
    logZ0tilde(1) = log(exp(-An[0]-An[1]+Cn[0]+zero) 
                            +exp( An[0]-An[1]-Cn[0]+zero));
   
    // Recurse through non-boundary values
  
    for (int n=3;n<r[0]+1;n++){       
        logZ1tilde(n-1) = log(exp(An[n-1]-Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(An[n-1]+Cn[n-2]+logZ1tilde(n-2)));
        logZ0tilde(n-1) = log(exp(-An[n-1]+Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(-An[n-1]-Cn[n-2]+logZ1tilde(n-2)));
    
    }
    
    // Compute marginal probability
    
    // Set to log(1/Z)
    logMargProb=-logZ;

    // Add log[\tilde{Z}_{r}(x_r)]
    if (x_r_rPLUSs[0]>0){ //X(r)=1
        logMargProb = logMargProb + logZ1tilde(r[0]-1);
    }else {              // X(r)=0
        logMargProb = logMargProb + logZ0tilde(r[0]-1);
    }

    // Add Z_{r+s}(x_{r+s})
    if (x_r_rPLUSs[s[0]]>0){ //X(r+s)=1
        logMargProb = logMargProb + logZ1(r[0]+s[0]-1);
    }else {               //X(r+s)=0
        logMargProb = logMargProb + logZ0(r[0]+s[0]-1);
    }

    // Add by \sum_{n=r}^{r+s-1} log[\phi(x_n,x_{n+1})]
    if (s[0]>0){ // otherwise this product is empty
        if (r[0]>1){ // no boundary \phi term needed
            for (int iter=1;iter<s[0]+1;iter++){ 
                n=iter+r[0]-1;
                logMargProb=logMargProb + (An[n]+Cn[n-1]*(2*x_r_rPLUSs[iter-1]-1))*(2*x_r_rPLUSs[iter]-1);
            }
                        
        }else{ // must compute boundary \phi term since r=1
            
            // n=1;
            logMargProb=logMargProb +  (2*x_r_rPLUSs[0]-1)*An[0]
                            +(2*x_r_rPLUSs[1]-1)*(An[1] +Cn[0]*(2*x_r_rPLUSs[0]-1));
            
            for (int iter2=2;iter2<s[0]+1;iter2++){ 
                logMargProb=logMargProb + (An[iter2]+Cn[iter2-1]*(2*x_r_rPLUSs[iter2-1]-1))*(2*x_r_rPLUSs[iter2]-1);
            } 
        }
    }
    
    // Write output to matlab
    logMargProbDoub[0] = logMargProb.toDouble();

}