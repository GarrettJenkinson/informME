// informME: An information-theoretic pipeline for WGBS data
// Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu), 
// and Jordi Abante (jabante1@jhu.edu)
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
 * Last Modified: 8/15/17
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
    
    // Make logZ0, logZ1 vectors 
    VectorXmp logZ1(N);
    VectorXmp logZ0(N);
    VectorXmp logZ0tilde(N);
    VectorXmp logZ1tilde(N);
   
    // Make mpreal scalars
    mpreal logZ;
    mpreal logMargProb;


    //
    // Load the cell arrays of strings and convert to mpreals
    // 
    mxArray *cell_element_ptr;

    char buf[2*(digits+2)]; // make a stack-allocated buffer with sufficient 
                            // space for string. theoretically digits+2 is 
                            // sufficient, but let's be paranoid to prevent troubles.

    for (int n=0;n<N;n++){
        //
        // get data for Z1
        //
        cell_element_ptr = mxGetCell(prhs[3],n); // cell_element_ptr is a pointer to the matlab string 
        // copy string into buffer buf
        if (mxGetString(cell_element_ptr, buf,(mwSize)(mxGetNumberOfElements(cell_element_ptr) + 1)) != 0)
            mexErrMsgTxt("Could not convert string data.");
        logZ = buf; // convert string to mpreal (done in mpreal package with equal sign)
        logZ1(n) = logZ; // store in vector

        //       
        // get data for Z0 
        //
        cell_element_ptr = mxGetCell(prhs[4],n); // cell_element_ptr is a pointer to the matlab string 
        // copy string into buffer buf
        if (mxGetString(cell_element_ptr, buf,(mwSize)(mxGetNumberOfElements(cell_element_ptr) + 1)) != 0)
            mexErrMsgTxt("Could not convert string data.");
        logZ = buf; // convert string to mpreal (done in mpreal package with equal sign)
        logZ0(n) = logZ; // store in vector

        //
        // get data for Z1tilde
        //
        cell_element_ptr = mxGetCell(prhs[6],n); // cell_element_ptr is a pointer to the matlab string 
        // copy string into buffer buf
        if (mxGetString(cell_element_ptr, buf,(mwSize)(mxGetNumberOfElements(cell_element_ptr) + 1)) != 0)
            mexErrMsgTxt("Could not convert string data.");
        logZ = buf; // convert string to mpreal (done in mpreal package with equal sign)
        logZ1tilde(n) = logZ; // store in vector

        //       
        // get data for Z0tilde
        //
        cell_element_ptr = mxGetCell(prhs[7],n); // cell_element_ptr is a pointer to the matlab string 
        // copy string into buffer buf
        if (mxGetString(cell_element_ptr, buf,(mwSize)(mxGetNumberOfElements(cell_element_ptr) + 1)) != 0)
            mexErrMsgTxt("Could not convert string data.");
        logZ = buf; // convert string to mpreal (done in mpreal package with equal sign)
        logZ0tilde(n) = logZ; // store in vector

    }

    cell_element_ptr = mxGetCell(prhs[5], 0); // cell_element_ptr is a pointer to the matlab string 
    // copy string into buffer buf
    if (mxGetString(cell_element_ptr, buf,(mwSize)(mxGetNumberOfElements(cell_element_ptr) + 1)) != 0)
        mexErrMsgTxt("Could not convert string data.");
    logZ = buf; // convert string to mpreal (done in mpreal package with equal sign)



    //    
    // Compute marginal probability
    //

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
