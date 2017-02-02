The following files are located in this directory:

FastaToCpGloc.m        - Script for analyzing the reference genome for 
                         CpG locations and for computing their densities 
                         and distances. This script should be run ONLY   
                         once before analyzing the .bam files. 

MatrixFromBAMfile.m    - Main function that takes a .bam file and creates 
                         all methylation information for a given chromosome. 
                         The output is written in separate directories for
                         each region in the form of methylation matrices.
						 WARNING: This function requires a working "samtools"
						 installation that is on the system $PATH.

nDensity.m             - Function that computes the density of the n-th 
                         CpG site. Called by FastaToCpGloc.m.

findSortedIndices.m    - Function that uses binary search to rapidly find 
                         the smallest and the largest index of the elements 
                         of a vector of sorted values within a a closed 
                         numerical region. Called by nDensity.m and 
                         MatrixFromBAMfile.m.

MatrixFromReads.m      - Function that takes SAM-formatted reads for a 
                         given region and constructs the matrix that 
                         stores the methylation information contained 
                         in the reads. Called by MatrixFromBAMfile.m. 

MergeMatrices          - Function to merge the results of parallel runs of 
                         MatrixFromBAMfile.m.


It also contains the following sub-directories:

genome                 - folder that stores the reference genome as well as
                         the results from FastaToCpGloc.m

matrices               - folder that stores the methylation matrices 
                         computed by MatrirxFromBAMfile.m
