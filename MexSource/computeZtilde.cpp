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
// computeZtilde.cpp

#include "mex.h"
#include "mpreal.h"
#include <iostream>
#include <Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	/* Recursive computation of the Ztilde function associated 
	 * with the 1D Ising model
     	 *
     	 * Author: W. Garrett Jenkinson
	 * 
     	 * Last Modified: November 18, 2016
     	 *
     	 * usage:
	 * [logZ1tilde,logZ0tilde,logZtilde]=computeZ(An,Cn);
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
    const double *An  = (double *)mxGetData(prhs[0]); // read in first rhs argument
    const double *Cn = (double *)mxGetData(prhs[1]); // read in second rhs argument
    
    plhs[0] = mxCreateCellMatrix(N, 1); // create first output argument in matlab
    plhs[1] = mxCreateCellMatrix(N, 1); // create second output argument in matlab
    plhs[2] = mxCreateCellMatrix(1, 1); // create third output argument in matlab
    
    // Make logZ0tilde and logZ1tilde vectors 
    VectorXmp logZ0tilde(N);
    VectorXmp logZ1tilde(N);
    
    // Make logZtilde scalar
    mpreal logZtilde;

    // Calculate first two boundary values
    logZ1tilde(0) = 0.0;
    logZ0tilde(0) = 0.0;
     
    logZ1tilde(1) = log(exp(-An[0]+An[1]-Cn[0]+zero) 
                            +exp(An[0]+An[1]+Cn[0]+zero));
    
    logZ0tilde(1) = log(exp(-An[0]-An[1]+Cn[0]+zero) 
                            +exp(An[0]-An[1]-Cn[0]+zero));
    
    // Recurse through non-boundary values
   
    for (int n=3;n<N+1;n++){    
        logZ1tilde(n-1) = log(exp(An[n-1]-Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(An[n-1]+Cn[n-2]+logZ1tilde(n-2))); 
        logZ0tilde(n-1) = log(exp(-An[n-1]+Cn[n-2]+logZ0tilde(n-2)) 
                                +exp(-An[n-1]-Cn[n-2]+logZ1tilde(n-2)));

    }

    // Compute log(Ztilde_N(0)+Ztilde_N(1))
    logZtilde = log( exp(logZ0tilde(N-1)) + exp(logZ1tilde(N-1)) );
    
    // Write output to Matlab
    for (int n = 0; n < N; n++) {
        mxSetCell(plhs[0],n,mxCreateString(strdup(logZ1tilde(n).toString().c_str())));
        mxSetCell(plhs[1],n,mxCreateString(strdup(logZ0tilde(n).toString().c_str())));

    }
    
    mxSetCell(plhs[2],0,mxCreateString(strdup(logZtilde.toString().c_str())));

}

