The following files are located in this directory:


compile.m
		        	    Matlab script to compile the mex C++ files in the directory.
						Must have configured a Matlab compiler, as well as
						a working mpfr library installation. This file should
						be modified to point to the relevant libraries on your
						system. The resulting mex binary files should be copied to
						the "Modeling" directory. This only needs to be done once.

calcMargProb.cpp        
						Calculates the marginal probability of a contiguous set 
						of CpG sites taking specific values.

computeMCtransProbs.cpp
					    This function computes the non-homogeneous markov chain transition
						probabilities as well as the initial conditions of this chain, which
						together are an alternative represntation of the Ising Model.

computeZ.cpp	        
						This function recurses through the Z computation on the Ising model.

computeZtilde.cpp       
						This function recurses through the Ztilde computation on the Ising model.

computeAveLogLikelihood.cpp          
					    This function computes the average log likelihood of a set of
						observations for a given Ising model.



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

Suppose we want to work inside a given directory such as the MexCode folder (which we will generically abreivate with ‘~’, which should be replaced by the full file path). Then we need to make a subdirectory to store all the libraries. We will call this mpfr:

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
    

Make updated compile.m to point towards the appropriate libraries/headers:	

The compile.m file should contain the following:

libpath = '~/mpfr/lib';
headpath2 = '~/mpfr/include';
headpath1 = '~/eigen';

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeZ.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeZtilde.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeAveLogLikelihood.cpp')


Run matlab and compile the MEX code:

	cd ~
	matlab –nojvm –nodisplay -nosplash –r “compile”

The output will be compiled mex code with extension “.mexa64”. These files should be moved into the Modeling directory.

DEBUGGING OLDER COMPILERS:

A problem encountered with older gcc compilers is that the linker doesn’t necessarily find linked libraries in unusual locations (such as the gmp and mpfr libraries we just installed), even when the binaries were compiled with the correct flags to tell the linker where to look. If this is the case, the following error will occur when calling the mex function inside matlab:

Invalid MEX-file '/cis/project/epigene/MexCode/computeZ.mexa64': libmpfr.so.4:
cannot open shared object file: No such file or directory

To correct this, add the following to your .bashrc file:

if [ -z "${LD_LIBRARY_PATH}" ]; then
    export LD_LIBRARY_PATH="~/mpfr/lib"
else
    export LD_LIBRARY_PATH="~/mpfr/lib:${LD_LIBRARY_PATH}"
fi

Generally, you must ensure that your gcc compiler is compatible with your version of matlab. Type "help mex" in matlab for more info.

