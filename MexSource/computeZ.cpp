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
// computeZ.cpp

#include "mex.h"
#include "mpreal.h"
#include <iostream>
#include <Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	/* Recursive computation of the partition function of 
	 * the 1D Ising model
	 *
	 * Author: W. Garrett Jenkinson
	 * 
     	 * Last Modified: November 18, 2016
     	 *
     	 * usage:
	 * [logZ1,logZ0,logZ]=computeZ(An,Cn);
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
    const double *An  = (double *)mxGetData(prhs[0]); // read in first rhs argument
    const double *Cn = (double *)mxGetData(prhs[1]); // read in second rhs argument
    
    plhs[0] = mxCreateCellMatrix(N, 1); // create first output argument in matlab
    plhs[1] = mxCreateCellMatrix(N, 1); // create second output argument in matlab
    plhs[2] = mxCreateCellMatrix(1, 1); // create third output argument in matlab
    
    // Make logZ0 and logZ1 vectors 
    VectorXmp logZ1(N);
    VectorXmp logZ0(N);
    
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
    
    // Compute log partition function = log(Z_1(0)+Z_1(1))
    logZ = log(exp(logZ0(0)) + exp(logZ1(0)));
    
    // Write output to Matlab

    for (int n = 0; n < N; n++) {
        mxSetCell(plhs[0],n,mxCreateString(strdup(logZ1(n).toString().c_str())));
        mxSetCell(plhs[1],n,mxCreateString(strdup(logZ0(n).toString().c_str())));

    }
    
    mxSetCell(plhs[2],0,mxCreateString(strdup(logZ.toString().c_str())));
    
}

