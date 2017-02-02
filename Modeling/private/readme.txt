This folder contains the following files:

EstimateParams.m 
            Main script that performs statistical estimation of parameters
            for the MRF model of methylation patterns (requires the Matlab 
            global optimization toolbox and symbolic math toolbox). 


procMatr.m
            function that preprocesses the data matrix to break apart
            non-contiguous reads

calcMargProb.m  
            function computes the marginal probabilities of a contiguous
            subset of the state vector 

computeAnCnm.m 
            function that computes the the a_n and C_{n,n+1} 
            parameters from their parameteric form

computeAveLogLikelihoo.m 
            function that computes the average log likelihood of the data
            matrix for a given model

computeMCtransProbs.m
            function that computes the non-homogenious markov chain 
            transition probabilities that can be used to represent the MRF

computeZ.m
            function that computes the partion function of the MRF using
            the forward recursion algorithm

computeZtilde.m
            function that computes the partition function of the MRF using
            the reverse recursion algorithm

objFnToMinimize.m
            function that computes the objective function that simulated
            annealing must minimize in order to find the Maximum Likelihood
            estimator of the \theta parameter vector

expectValSXexact.m
            function to compute the expected value of S(X) 

runEstimation.sh
            Shell script to submit a single job to Sun Grid Engine cluster

runAllInitEstimation.sh
            Shell script to submit multiple jobs to a Sun Grid Emgine 
            using "qsub" commands

runAllEstimation.m
            matlab script to run estimation code

matrixNumbers.txt
            file containing the numbers of the matrices that needs to be
            fit initially in the first round of clustering. File is used
            by runAllInitEstimation.sh


It also contains the following subdirectory:

results   -- a folder into which EstimateParams.m writes its results

optiAlgs -- a folder containing code
            [modified to ensure compatability with modern Matlab software] of 
            the MCS optimization algorithm as well as the  MINQ (bound constrained 
            quadratic program solver) and GLS (global line search) software that this 
            package depends on. The original software of these three packages along
            with their corresponding documentation as downloaded from 
            http://www.mat.univie.ac.at/~neum/software/mcs/
            http://www.mat.univie.ac.at/~neum/software/minq/
            http://www.mat.univie.ac.at/~neum/software/ls/
            is stored in the "OptimizationAlgorithms" subdirectory. 