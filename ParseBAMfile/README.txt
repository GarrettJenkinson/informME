Last Modified: 12/10/16

----------------------------------------------------------------
DIRECTORY ParseBAMfile
----------------------------------------------------------------

----------------------------------------------------------------
Scripts to run these files on a Sun Grid Engine cluster can be 
found in ../SubmissionScripts/ParseBAMfile. 
----------------------------------------------------------------

This directory includes the following files:

FastaToCpG.m
------------

This function is used to analyze a reference genome in order to 
find and store the locations of all CpG sites within each 
chromosome and to compute the CpG densities at each CpG site 
as well as the distances between neighboring CpG sites. A 1-based 
coordinate system is used, in which the first base is assigned 
to position 1 and the location of a CpG site is defined by the 
position of the C nucleotide on the forward strand of the 
reference genome. 

findSortedIndices.m
-------------------

This function uses binary search to rapidly find the smallest 
and the largest indices of the elements of a vector of sorted 
values (low to high) that are between a lower and an upper 
bound.

MatrixFromBAMfile.m
-------------------

This function processes a BAM file with aligned reads to a 
reference genome and produces methylation information for 
nonoverlapping genomic regions (containing the same number 
of base pairs) in a given chromosome. The final output for 
each genomic region is a matrix with -1,0,1 values. Each row 
of the matrix is a methylation read, whereas each column 
represents a CpG site within the genomic region. A value 
of -1 indicates no methylation information is available for 
the CPG site, 0 indicates that the CpG site is unmethylated, 
and 1 indicates that the CpG site is methylated. 

MatrixFromReads.m
-----------------

This function determines the methylation status of the CpG 
sites within a given genomic region using SAM-formatted reads 
with base pair coverage within the region. The output is a 
matrix with -1,0,1 values. Each row of the matrix is a single 
sequencing read, whereas each column represents a CpG site 
within the genomic region. A value of -1 indicates no 
methylation information is available for the CPG site, 
0 indicates that the CpG site is unmethylated, and 1 
indicates that the CpG site is methylated. 

MergeMatrices
-------------

This function merges the outputs from MatrixFromBAMfile.m 
into a single MATLAB MAT file comprised of a single hashtable.

nDensity.m
----------

This function computes the density of the n-th CpG site 
within a given chromosome. 
              

This directory also contains the following subdirectories:

genome
------               

This subdirectory contains the reference genome FASTA file 
as well the results obtained by FastaToCpG.m. The results 
are organized with a subdirectory for each species (Human, 
Mouse, etc.), with each species folder including a MATLAB 
MAT file CpGlocationChr#.mat for each chromosome that contains 
the following information: 
o location of CpG sites 
o CpG density for each CpG site 
o distance between neighboring CpG sites
o location of the last CpG site in the chormosome
o length of chromosome (in # of base pairs)

indexedBAMfiles
---------------

This subdirectory contains indexed BAM files (.bam and .bai) 
with reads aligned to a reference genome that are to be 
processed in order to produce methylation information about 
the status of CpG sites in the genome. These files are 
organized within a subdirectory for each species (Human, 
Mouse, etc.).

matrices
--------

This subdirectory contains methylation matrices, obtained 
by processing available BAM files, organized with a 
subdirectory for each species (Human, Mouse, etc.), with 
each species folder including subdirectories for each 
chromosome (chr1, chr2, etc.). Each chromosome subdirectory 
includes the MATLAB file BAMfileName.mat for each BAM 
file processed, with each file containing the following 
information for each genomic region to be subsequently 
used in model estimation: 
o data matrix with -1,0,1 values for methylation status
o CpG locations broken down by region.

