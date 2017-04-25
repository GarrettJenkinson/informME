Last Modified: 04/22/17

----------------------------------------------------------------
DIRECTORY MexSource
----------------------------------------------------------------

----------------------------------------------------------------
Scripts to run these files on a Sun Grid Engine cluster can be 
found in ../SubmissionScripts/MexSource. 
----------------------------------------------------------------

This directory includes the following files:

calcMargProb.cpp
----------------

This function computes the joint probability of a 
contiguous set of CpG sites using Equation (25) in [1].

compile.m
---------					

This function creates MEXfiles from the C++ files 
located in the directory "MexSource". This requires 
configuration of a MATLAB compiler, as well as 
working installations of the MPFR (http://www.mpfr.org) 
and EIGEN (eigen.tuxfamily.org) packages. The generated 
binary MEX files are automatically placed in the 
"Modeling" and "SingleAnalysis" directories of informME.

computeAveLogLikelihood.cpp
---------------------------

This function computes the average "marginalized" 
log-likelihood function of a set of methylation 
observations for a given Ising model, given by 
Equation (27) in [1].

computeMCtransProbs.cpp
-----------------------

This function computes the initial and transition 
probabilities of the alternative inhomogeneous Markov 
chain representation of the 1D Ising model given by 
Equations (17) and (18) in [1].

computeZ.cpp
------------
	        
This function recursively computes the partition 
function of the 1D Ising model using Equation (15) 
in [1].

computeZtilde.cpp 
-----------------
      
This function recursively computes the Ztilde 
function associated with the 1D Ising model using 
Equation (24) in [1].


REFERENCES
----------

[1] Jenkinson, G., Feinberg, A.P., and Goutsias, J. (2017) 
    An information-theoretic approach to the modeling and    
    analysis of whole-genome bisulfite sequencing data, 
    Submitted.


DEPENDENCIES
------------

The C++ code has the following dependencies:

eigen: The Eigen C++ library for linear algebra. 
       Tested with v3.2.1.

GMP: The GNU Multiprecision library. 
     Tested with v5.1.3.

MPFR: The GNU MPFR library for multiprecision 
      arithmetic. Tested with v.3.1.2.

MPFR C++: The MPFR C++ header library by Pavel 
Holoborodko. Tested with v.3.5.7.
