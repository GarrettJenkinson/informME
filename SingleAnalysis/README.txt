Last Modified: 02/03/17

----------------------------------------------------------------
DIRECTORY SingleAnalysis
----------------------------------------------------------------

----------------------------------------------------------------
Scripts to run these files on a Sun Grid Engine cluster can be 
found in ../SubmissionScripts/SingleAnalysis. 
----------------------------------------------------------------

This directory includes the following files:

computeAnCn.m
-------------

This function computes the a_n and c_n parameters of the Ising 
model within a genomic region by using the CpG densities, the 
CpG distances, and the estimated alpha, beta and gamma 
parameters of the model. 

computeLstats.m
---------------

This function computes the probability distribution of the 
methylation level L within a genomic region comprised N CpG 
sites as well as the corresponding normalized methylation 
entropy.

exactSampling.m
---------------

This function draws an exact Monte Carlo sample from the 
Ising model within a genomic region containing N CpG sites.

findSortedIndices.m
-------------------

This function uses binary search to rapidly find the smallest 
and the largest indices of the elements of a vector of sorted 
values (low to high) that are between a lower and an upper 
bound.

h_func.m
--------

This function computes the (log2-based) entropies of a 
collection X(q,r), q = 1,2,...,Q, r = 1,2,...,R of QxR binary 
random variables with corresponding probabilities p(q,r) 
and 1-p(q,r). 

MakeBEDsForMethAnalysis.m
-------------------------

This function makes BED files for the methylation analysis 
results obtained by means of MethAnalysisForChr.m for a 
single phenotype. 

MergeMethAnalysis.m
-------------------

This function merges the output of MethAnalysisForChr.m into 
a single hashtable.

MethAnalysisForChr.m
---------------------

This function performs methylation analysis of a given 
chromosome in a single phenotype. The function can be used 
on a computing cluster to break the analysis work to many 
independent parallel job processes. This is performed only 
after EstParamsForChr.m in the Modeling subdirectory is run 
to build the Ising models for the phenotype.

MethAnalysisForRegion.m
-----------------------

This function performs methylation analysis of a genomic 
region used for estimating the parameters of the Ising 
model by computing a number of statistical summaries of 
the methylation state within the region, including 
probability distributions of methylation levels, mean 
meathylation levels, and normalized methylation entropies. 
If desired, this function also computes entropic sensitivity 
indices, as well information-theoretic quantities associated 
with methylation channels, such as turnover ratios, channel 
capacities, and relative dissipated energies. 


This directory also contains the following subdirectories:

private
-------

This subdirectory contains the MATLAB code that implements 
the maxEnt method presented in [1]. The code has been 
modified to ensure compatability with new versions of 
MATLAB. This is a non-GPL code, whose author has provided 
permission to modify and distribute with the informME 
package. Parties interested in using this code outside 
of informME should contact the author of the code directly.

results
-------

This subdirectory contains the single-sample methylation 
analysis results organized within a subdirectory for 
each species (Human, Mouse, etc.), with each species 
folder containing a subdirectory for each chromosome 
(chr1, chr2, etc.). Each chromosome subdirectory 
contains the MATLAB file phenoName_Analysis.mat for 
each phenotypic methylation sample (colonnormal, 
coloncancer, etc.), with each phenoName_Analysis.mat 
file containing the following information for each 
genomic region used in model estimation: 
o The locations of the CpG sites within the genomic region.
o Numbers of CpG sites within the analysis subregions. 
o Which analysis subregions are modeled and which are not.
o Estimated parameters of Ising model in genomic region
o Methylation level probabilities in modeled subregions.
o Coarse methylation level probabilities.
o Mean methylation levels.
o Normalized methylation entropies.
o Entropic sensitivity indices (if ESIflag = 1).
o Turnover ratios (if MCflag = 1).
o Channel capacities (if MCflag = 1).
o Relative dissipated energies (if MCflag = 1).


OUTPUT
------

SingleAnalysis generates the following methylation tracks 
of single-sample analysis in the form of BED files, which 
are written within the subdirectory /informME/BEDfiles 
organized for each species (Human, Mouse, etc.):

o MML-phenoName.bed
      mean methylation levels

o NME-phenoName.bed
      normalized methylation entropy

o METH-phenoName.bed
       methylation-based classification (non-variable)

o VAR-phenoName.bed
      methylation-based classification (variable)

o ENTR-phenoName.bed
       entropy-based classification

o ESI-phenoName.bed (if ESIflag = 1)
      entropic sensitivity indices

o TURN-phenoName.bed (if MCflag = 1) 
       turnover ratios

o CAP-phenoName.bed (if MCflag = 1) 
      channel capacities

o RDE-phenoName.bed (if MCflag = 1) 
      relative dissipated energies


REFERENCES
----------

[1] Mohammad-Djafari, A. (1992) A Matlab program to calculate 
    the maximum entropy distributions. In Fundamental Theories 
    of Physics: Maximum Entropy and Bayesian Methods, C.R. Smith 
    et al. (eds.), 50:221-233, Kluwer Academic Publishers.
