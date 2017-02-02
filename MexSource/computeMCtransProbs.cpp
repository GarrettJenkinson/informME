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
#include <sstream>
#include <iostream>
#include <Eigen/MPRealSupport>
#include "mpreal.h"
#include "mex.h"

using namespace mpfr;
using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/* Recursive Computation of Markov Chain Transition Probabilities
     *
     * Authors: W. Garrett Jenkinson
	 * 
     * Last Modified: March 2, 2014
     *
     * usage:
	 * [p1,transProbs] = computeMCtransProbs(invOmega_An,invOmega_Cnm,logZ1,logZ0,logZ);
     *
     INPUTS:
     Omega_alpha is a Nx1 vector of Omega*alpha_n parameters for n=1,...,N
     z1 is a Nx1 vector of values Z_n(x_n) for x_n=1 and n=1,...,N
     z0 is a Nx1 vector of values Z_n(x_n) for x_n=0 and n=1,...,N
     Z is the partition function
    
     OUTPUTS:
     p1 is the marginal probability  P(X_1=0)
     transProbs is a (N-1)x2 matrix
       -First  column has probability P(X_{n+1}=0|X_n=0) for n=1,2,...,N-1
       -Second column has probability P(X_{n+1}=0|X_n=1) for n=1,2,...,N-1
    
     Note the joint dist can be calculated from these probabilities:
     P(\bfX=\bfx) = P(X_1)\prod_{n=1}^{N-1} P(X_{n+1|X_n)
    
     More importantly these probabilities can be used to iteratively draw an
     exact sample from the joint distribution without MCMC.
     */
	 
    
    /*
    // Setup Multiple Precision Code
    */
    
    using mpfr::mpreal;
    
    // Required precision of computations in decimal digits
    // Play with it to check different precisions
    const int digits = 200;
    
    // Initialize Multiple Precision Number 
    const mpreal zero = 0.0;
    
    // Setup default precision for all subsequent computations
    // MPFR accepts precision in bits - so we do the conversion
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    
    // Declare matrix and vector types with multi-precision scalar type
    typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
    typedef Matrix<mpreal,Dynamic,1>        VectorXmp;
    
    
    /*
	//Initialize matlab input/output vars
    */
	
    int N = mxGetNumberOfElements(prhs[0]); //number of elements of first rhs argument
    const double *invOmega_An  = (double *)mxGetData(prhs[0]);  // read in first rhs argument
    const double *invOmega_Cnm = (double *)mxGetData(prhs[1]);  // read in second rhs argument
    
    //mwSize ndims[] = {N-1,2};
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N-1,2,mxREAL);
    //plhs[0] = mxCreateCellMatrix(1, 1);  // create first output argument in matlab
    //plhs[1] = mxCreateCellArray(2, ndims);// create second output argument in matlab
    
    //Make mpreal outputs (later converted to strings when passed to matlab) 

    MatrixXmp transProbs(N-1,2);
   
    mpreal p1;

    //Make logZ0 and logZ1 vectors 

    VectorXmp logZ0(N);
    VectorXmp logZ1(N);
    
    //Make logZ scalar
    mpreal logZ;
    
    /*
    // Read in logZ0, logZ1, logZ from matlab strings and cast as mpreal values
    */

    /*
    mxArray *cell_element_ptr;
    string   buf;
    int      buflen,status;
    
    for (int n = 0; n < N; n++) {
        // first do logZ1
        cell_element_ptr = mxGetCell(prhs[2], n); // get pointer to cell string
        //buflen = (mxGetM(cell_element_ptr) * mxGetN(cell_element_ptr)) + 1;//Find out how long the input string array is. 
        buf = mxArrayToString(cell_element_ptr);
        logZ1(n)=buf;   
        
        
        // next do logZ0
        cell_element_ptr = mxGetCell(prhs[3], n); // get pointer to cell string
        //buflen = (mxGetM(cell_element_ptr) * mxGetN(cell_element_ptr)) + 1;//Find out how long the input string array is.  
        buf = mxArrayToString(cell_element_ptr);
        logZ0(n)=buf;    
    }
    
    // finally do logZ
    cell_element_ptr = mxGetCell(prhs[4], n); // get pointer to cell string
    //buflen = (mxGetM(cell_element_ptr) * mxGetN(cell_element_ptr)) + 1;//Find out how long the input string array is.   
    buf = mxArrayToString(cell_element_ptr);
    logZ1(n)=buf;   
    */
    
    
    /* 
     *Begin computation of logZ's
     */
    
    
    //Calculate first boundary values
    logZ1(N-1)=0.0;
    logZ0(N-1)=0.0;
    
    //Recurse through non-boundary values
    for (int n=N-1;n>1;n--){
        //log(Z_n(x_n)) = log( exp( log(\phi_n(x_n,0))+log(Z_{n+1}(0) ) ...
        //                     +exp( log(\phi_n(x_n,1))+log(Z_{n+1}(1) )) ;
        logZ1(n-1) = log(  exp(-invOmega_An[n]-invOmega_Cnm[n-1]+logZ0(n))
                         + exp( invOmega_An[n]+invOmega_Cnm[n-1]+logZ1(n))  );
        
        
        
        logZ0(n-1) = log( exp(-invOmega_An[n]+invOmega_Cnm[n-1]+logZ0(n))
                        + exp( invOmega_An[n]-invOmega_Cnm[n-1]+logZ1(n))  );
        
        
    }
    
    //Calculate last boundary values
    
    logZ1(0) = log( exp( invOmega_An[0]-invOmega_An[1]-invOmega_Cnm[0]+logZ0(1))
                  + exp( invOmega_An[0]+invOmega_An[1]+invOmega_Cnm[0]+logZ1(1)) );
    
    
    
    logZ0(0) = log( exp(-invOmega_An[0]-invOmega_An[1]+invOmega_Cnm[0]+logZ0(1))
                  + exp(-invOmega_An[0]+invOmega_An[1]-invOmega_Cnm[0]+logZ1(1)) );
    
    
    //Compute log partition function = log(Z_1(0)+Z_1(1))
    logZ = log( exp(logZ0(0)) + exp(logZ1(0)) );
    
    
    

    /*
    % Compute boundary condition 
    */
    
    p1 = exp( logZ0(0)-logZ ); // P(X_N=0)
    

    /*
    % P(x_2=0|x_1=0) = \phi_1(x_1=0,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=0)
    */
    
    transProbs(0,0) = exp( -invOmega_An[0]-invOmega_An[1]+invOmega_Cnm[0]+logZ0(1)-logZ0(0) ); 

            
    /*
     * P(x_2=0|x_1=1) = \phi_1(x_1=1,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=1)
     */
    
    transProbs(0,1) = exp( invOmega_An[0]-invOmega_An[1]-invOmega_Cnm[0]+logZ0(1)-logZ1(0) );
     
    // debugging
    /*
    stringstream s;
    s << "transProbs(0,1)=" << transProbs(0,1).toDouble() << "\n";
    string stri; 
    stri = s.str();
    mexPrintf(stri.c_str());
    */

    /*
     * Use backwards recursion to compute transition probabilities
     */
    for (int n=2;n<N;n++){
        /*
        % P(x_{n+1}=0|x_n=0) = \phi_n(x_n=0,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=0)
        */
        
        transProbs(n-1,0) = exp( -invOmega_An[n]+invOmega_Cnm[n-1]+logZ0(n)-logZ0(n-1) );
        
        

        /*
        % P(x_{n+1}=0|x_{n}=1) = \phi_n(x_n=1,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=1)
        */
         
        transProbs(n-1,1) = exp( -invOmega_An[n]-invOmega_Cnm[n-1]+logZ0(n)-logZ1(n-1) );
    }

    
    /*
     *Write output to Matlab
     */
    
    
    /*
    mxSetCell(plhs[0],0,mxCreateString(strdup(p1.toString().c_str())));
    int nsubs=2;int arrIndex;
    mwIndex subs[2];
    for (int n = 0; n < N-1; n++) {
        //set first column of trans probs
        subs[0]=n;subs[1]=0;
        arrIndex = mxCalcSingleSubscript(plhs[1], nsubs, subs);
        mxSetCell(plhs[1],arrIndex,mxCreateString(strdup(transProbs(n,0).toString().c_str()))); 
        
        //set second column of trans probs
        subs[1]=1;
        arrIndex = mxCalcSingleSubscript(plhs[1], nsubs, subs);
        mxSetCell(plhs[1],arrIndex,mxCreateString(strdup(transProbs(n,1).toString().c_str()))); 
    }*/

    double *outputProb;
    outputProb = mxGetPr(plhs[0]);
    outputProb[0]=p1.toDouble();
    
    
    double *outputMatrix;
    outputMatrix = mxGetPr(plhs[1]);
    for (int n = 0; n < N-1; n++) {
        outputMatrix[n]   = transProbs(n,0).toDouble();
        outputMatrix[n+N-1] = transProbs(n,1).toDouble();
    }
    

}
