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
    analysis of whole-genome bisulfite sequencing data, 
    Bioinformatics, XX, XXX-XXX.

The C++ code has the following dependencies:

eigen                   
						The Eigen C++ library for linear algebra. Tested with v3.2.1.

GMP
                        The GNU Multiprecision library. Tested with v5.1.3.

MPFR
                        The GNU MPFR library for multiprecision arithmetic. Tested with v.3.1.2.

MPFR C++
                        The MPFR C++ header library by Pavel Holoborodko. Tested with v.3.5.7.


------------------------------------------------------------------------------------------------
Here are general instructions for installing the dependencies 
(intended as a guide, links may break, details might change slightly, etc.):
------------------------------------------------------------------------------------------------
Install GMP and MPFR libraries:

Suppose we want to work inside a given directory such as the MexCode folder (which we will 
generically abreivate with ‘~’, which should be replaced by the full file path). Then we 
need to make a subdirectory to store all the libraries. We will call this mpfr:

	mkdir mpfr

Then we create a directory tree similar to be found in /usr/ 
	
	cd mpfr
	mkdir bin
	mkdir include
	mkdir lib
	mkdir lib32
	mkdir lib64
	mkdir local
	mkdir sbin
	mkdir share
	mkdir src

Then we download the latest versions of the GMP and MPFR in the src subdirectory:

	cd src
	wget http://mirror.team-cymru.org/gnu/gmp/gmp-5.1.3.tar.gz
	wget http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz

Next we expand these archive files:

	tar xzf gmp-5.1.3.tar.gz
	tar xzf mpfr-3.1.2.tar.gz

Now we go into the gmp folder and install gmp from source:

	cd gmp-5.1.3
	./configure --prefix=~/mpfr
	make
	make check
	make install

Next we will go into the mpfr folder and install mpfr from source:
	
	cd ../mpfr-3.1.2
	./configure --prefix=~/mpfr --with-gmp=~/mpfr	
	make
	make check
	make install


Download the Eigen Header Library:

	cd ~
	wget --no-check-certificate http://bitbucket.org/eigen/eigen/get/3.2.1.tar.gz
	tar xzf 3.2.1.tar.gz
	mv eigen-eigen-6b38706d90a9 eigen

Copy the “unsupported” MPRealSupport to the Eigen folder:

	cd eigen
	cp unsupported/Eigen/MPRealSupport Eigen/MPRealSupport		

Download mpreal.h from the author of MPFR C++ and place it in the eigen directory.

    cd ~
    wget --no-check-certificate http://www.holoborodko.com/pavel/wp-content/plugins/download-monitor/download.php?id=4
    unzip mpfrc*.zip -d ~/mpfrc
    mv ~/mpfrc/mpreal.h ~/eigen
    

Now just run compile.m in matlab per the instructions in that file to compile the C++ code.

DEBUGGING OLDER COMPILERS:

A problem encountered with older gcc compilers is that the linker doesn’t necessarily find linked libraries 
in unusual locations (such as the gmp and mpfr libraries we just installed), even when the binaries were 
compiled with the correct flags to tell the linker where to look. If this is the case, the following error 
will occur when calling the mex function inside matlab:

Invalid MEX-file '~/MexCode/computeZ.mexa64': libmpfr.so.4:
cannot open shared object file: No such file or directory

To correct this, add the following to your .bashrc file:

if [ -z "${LD_LIBRARY_PATH}" ]; then
    export LD_LIBRARY_PATH="~/mpfr/lib"
else
    export LD_LIBRARY_PATH="~/mpfr/lib:${LD_LIBRARY_PATH}"
fi

Generally, you must ensure that your gcc compiler is compatible with your version of matlab. 
Type "help mex" in matlab for more info.


