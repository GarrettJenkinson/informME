Last Modified: 12/03/16

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
contiguous set of CpG sites using Equations (S18) 
and (S19) in the Supplementary Methods of [1].

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
Equation (S26) in the Supplementary Methods of [1].

computeMCtransProbs.cpp
-----------------------

This function computes the initial and transition 
probabilities of the alternative inhomogeneous Markov 
chain representation of the 1D Ising model given by 
Equations (S15) and (S16) in the Supplementary Methods 
of [1].

computeZ.cpp
------------
	        
This function recursively computes the partition 
function of the 1D Ising model using Equation (S13) 
in the Supplementary Methods of [1].

computeZtilde.cpp 
-----------------
      
This function recursively computes the Ztilde 
function associated with the 1D Ising model using 
Equations (S22) and (S23) in the Supplementary 
Methods of [1].


REFERENCES
----------

[1] Jenkinson, G., Feinberg, A.P., and Goutsias, J. (2017) 
    informME: An information-theoretic pipeline for methylation 
    analysis of whole-genome bisulfite sequencing data, Submitted.


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
