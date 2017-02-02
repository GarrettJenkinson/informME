This folder contains the following files:

EstimateParams.m 
            Main function that performs statistical estimation of parameters
            for the MRF model of methylation patterns (requires the Matlab 
            symbolic math toolbox). 

AnalysisRegion.m
			This function computes the statistical summaries in a specified 
			genomic region.

MakeBEDfilesAnalysisComp.m
			This function makes the relevant bed files to compare two phenotypes. It
			assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
 			with AnalysisForChr.m/MergeAnalysis.m for both phenotypes.

MakeBEDfilesAnalysis.m
			This function makes the relevant bed files for a given phenotype. It
			assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
			with AnalysisForChr.m/MergeAnalysis.m for the phenotype.

EstParamsForChr.m
			This function takes a list of bam files (all corresponding to the 
			same phenotype) and performs statistical model estimation for regions 
			within the specified chromosome of interest. In particular, this 
			function can be used on a computing cluster to break up the work to 
			many independent parallel job processes. This is performed only after 
			MatrixFromBamfile.m has been run to produce the data formats required 
			for statistical esimtation.

AnalysisForChr.m
			This function takes a list of bam files (all corresponding to 
			the same phenotype) and performs statistical model estimation 
			for regions within the specified chromosome of interest. In 
			particular, this function can be used on a computing cluster to 
			break up the work to many independent parallel job processes. 
			This is performed only after MatrixFromBamfile.m has been run 
			to produce the data formats required for statistical esimtation.

calcMargProb.m  
            function computes the marginal probabilities of a contiguous
            subset of the state vector. Should use C++ mex compiled version 
			instead for speed. 

computeAnCnm.m 
            function that computes the the a_n and C_{n,n+1} 
            parameters from their parameteric form

computeAveLogLikelihood.m 
            function that computes the average log likelihood of the data
            matrix for a given model. Should use C++ mex compiled version 
			instead for speed. 

computeMCtransProbs.m
            function that computes the non-homogenious markov chain 
            transition probabilities that can be used to represent the MRF
			Should use C++ mex compiled version instead for speed. 

computeZ.m
            function that computes the partion function of the MRF using
            the forward recursion algorithm. Should use C++ mex compiled  
			version instead for speed. 

computeZtilde.m
            function that computes the partition function of the MRF using
            the reverse recursion algorithm. Should use C++ mex compiled  
			version instead for speed. 

objFnToMinimize.m
            function that computes the objective function that simulated
            annealing must minimize in order to find the Maximum Likelihood
            estimator of the \theta parameter vector

expectValSXexact.m
            function to compute the expected value of S(X) 

calcMargProb.m
			Calculates the marginal probability of a contiguous set of CpG 
			sites taking specific values. Should use C++ mex compiled  
			version instead for speed. 

computeYprobs.m
			This function computes the probability distribution of level
			Y=(1/N)*\sum_n X_n

DiffEntrBlocks.m
			This function does postprocessing of bed files to find DEBs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were 
			analyzed with AnalysisForChr.m/MergeAnalysis.m for both 
			phenotypes and bed files generated with MakeBEDfilesAnalysisComp.m.

DiffEntrRegions.m
			This function does postprocessing of bed files to find DERs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were 
			analyzed with AnalysisForChr.m/MergeAnalysis.m for both 
			phenotypes and bed files generated with MakeBEDfilesAnalysisComp.m.

DiffMethBlocks.m
			This function does postprocessing of bed files to find DMBs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were 
			analyzed with AnalysisForChr.m/MergeAnalysis.m for both 
			phenotypes and bed files generated with MakeBEDfilesAnalysisComp.m.

DiffMethRegions.m			
			This function does postprocessing of bed files to find DMRs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were 
			analyzed with AnalysisForChr.m/MergeAnalysis.m for both 
			phenotypes and bed files generated with MakeBEDfilesAnalysisComp.m.

exactSampling.m
			A function that draws an exact Monte Carlo sample from the Ising model.

expectValSXexact.m
			This function computes the E[S(X)] as well as the nearest neighbor 
			correlation coefficients for the Ising model.

findSortedIndices.m
			This function uses binary search to rapidly find the smallest and the 
			largest index of the elements of a vector x of sorted values (from 
	 		low to high) within a closed numerical region.

h_func.m
			This function returns the entropy for input probabilities 0<=p<=1
		 	Note: for entries of p outside this range, h_func returns 0.

maxent.m
			File from 1992 paper "A Matlab program to calculate the maximum  
			entropy distributions" by Ali Mohammad-Djafari. Modified slighty by 
			Garrett Jenkinson with comments indicating the modifications.

MergeAnalysis.m
			This function merges the output of AnalysisForChr.m into a single
			hashtable.

MergeEstParams.m
			This function merges the output of EstParamsForChr.m into a single
			hashtable.

MethBlocks.m
			This function does PostProcessing of bed files to find MBs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
			with AnalysisForChr.m/MergeAnalysis.m  and bed files generated with
			MakeBEDfilesAnalysis.m.

MethRegions.m
			This function does PostProcessing of bed files to find MRs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
			with AnalysisForChr.m/MergeAnalysis.m  and bed files generated with
			MakeBEDfilesAnalysis.m.

OrdDisordBlocks.m
			This function does PostProcessing of bed files to find EBs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
			with AnalysisForChr.m/MergeAnalysis.m  and bed files generated with
			MakeBEDfilesAnalysis.m.

OrdDisordRegions.m
			This function does PostProcessing of bed files to find ERs. 
			It assumes that statistical estimation has been completed using
			EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
			with AnalysisForChr.m/MergeAnalysis.m  and bed files generated with
			MakeBEDfilesAnalysis.m.

processMatrix.m
			This function takes a data matrix and breaks apart non-contiguous
			observations of CpG sites into independent reads. It also provides the
			start and end index of each read, allowing the data matrix to be parsed
			more quickly in downstream applications. 

readme.txt
			See readme.txt for an infinite recursion.



It also contains the following subdirectory:

results   -- a folder into which EstimateParams.m writes its results

private -- a folder containing code licensed by their respective authors to Garrett
            Jenkinson to be modified and distributed with the informME software package.
            See the headers of these files for detailed licensing information.
            Thanks to the authors Arnold Neumaier and Ali Mohammad-Djafari for their
            kind licensing permissions. Contact the authors directly for licensing
            questions beyond usage with this package. 
