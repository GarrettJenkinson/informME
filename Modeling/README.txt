Last Modified: 12/08/16

----------------------------------------------------------------
DIRECTORY Modeling
----------------------------------------------------------------

----------------------------------------------------------------
Scripts to run these files on a Sun Grid Engine cluster can be 
found in ../SubmissionScripts/Modeling. 
----------------------------------------------------------------

This directory includes the following files:

computeAnCn.m
--------------

This function computes the a_n and c_n parameters of the Ising 
model within a genomic region by using the CpG densities, the 
CpG distances, and the estimated alpha, beta and gamma 
parameters of the model. 

EstimateParams.m
----------------

This function takes a list of BAM files (which correspond to the 
same phenotype) and estimates the parameters of the 1D Ising model 
that best fits the methylation data associated with a specific 
genomic region.

EstParamsForChr.m
-----------------

This function takes a list of BAM files (which correspond to the 
same phenotype) and performs statistical model estimation within 
a specific chromosome of interest. The function can be used on a 
computing cluster to break the work of model estimation to many 
independent parallel job processes. This is performed only after 
MatrixFromBAMfile.m in the ParseBAMfile subdirectory is run to 
produce the data required for statistical estimation.

findSortedIndices.m
-------------------

This function uses binary search to rapidly find the smallest 
and the largest indices of the elements of a vector of sorted 
values (low to high) that are between a lower and an upper 
bound.

MergeEstParams.m
----------------
			
This function merges the output of EstParamsForChr.m into a 
single hashtable.

objFnToMinimize.m
-----------------

This function evaluates the objective function that must be 
minimized in order to obtain maximum-likelihood estimates 
for the alpha, beta, and gamma parameters of the Ising model. 

processMatrix.m
---------------

This function takes a data matrix and breaks non-contiguous 
observations of CpG sites into independent reads. It also 
provides the start and end indices of each read, allowing 
the data matrix to be parsed more quickly in downstream 
applications.

	
This directory also contains the following subdirectories:

private
-------

This subdirectory has code which has license provided by
its author, Dr. Arnold Neumaier, to Garrett Jenkinson to 
distribute and modify as part of the informME software package. 
See the header of these files for more licensing information.
This subdirectory contains the MATLAB code that implements 
the Multivel Coordinated Search (MCS) optimization algorihm 
presented in [1]. It also contains the MINQ, a bound 
constrained quadratic program solver, and the global line 
search (GLS) software that MCS depends on. This code has 
been modified to ensure compatability with new versions of 
MATLAB. The original software of these three packages along 
with their corresponding documentation is available at:
http://www.mat.univie.ac.at/~neum/software/mcs/
http://www.mat.univie.ac.at/~neum/software/minq/
http://www.mat.univie.ac.at/~neum/software/ls/

results
-------

This subdirectory contains the parameter estimation 
(modeling) results organized with a subdirectory for each 
species (Human, Mouse, etc.), with each species folder 
containing a subdirectory for each chromosome (chr1, chr2, 
etc.). Each chromosome subdirectory contains the MATLAB 
file phenoName.mat for each phenotypic methylation sample 
(colonnormal, coloncancer, etc.), with each phenoName.mat 
file containing the following information for each genomic 
region used in model estimation: 
o CpG distances
o CpG densities
o estimated alpha, beta, and gamma parameters of the Ising 
  model
o initial and transition probabilities of the inhomogeneous 
  Markov chain representation of the Ising model
o marginal probabilities at each CpG site
o the log partition function of the estimated Ising model.

REFERENCES
----------

[1] Huyer, W. and Neumaier, A. (1999) Global optimization
    by multilevel coordinated search. J. Global Optim., 
    14:331-355.
