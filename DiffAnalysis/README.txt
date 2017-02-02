Last Modified: 12/09/16

----------------------------------------------------------------
DIRECTORY DiffAnalysis
----------------------------------------------------------------

----------------------------------------------------------------
Scripts to run these files on a Sun Grid Engine cluster can be 
found in ../SubmissionScripts/DiffAnalysis. 
----------------------------------------------------------------

This directory includes the following files:

findSortedIndices.m
-------------------

This function uses binary search to rapidly find the smallest 
and the largest indices of the elements of a vector of sorted 
values (low to high) that are between a lower and an upper 
bound.

MakeBEDsForDiffMethAnalysis.m
-----------------------------

This function makes BED files for the differential version 
of the methylation analysis results obtained by means of 
MethAnalysisForChr.m applied on two dinstict phenotypes. 


OUTPUT
------

The code in this directory generates the following tracks 
of differential methylation analysis in the form of BED 
files, which are written in the subdirectory /informME/BEDfiles 
organized for each species (Human, Mouse, etc.):

o dMML-phenoName1-VS-phenoName2.bed
       differences in mean methylation levels

o DMU-phenoName1-VS-phenoName2.bed
      differential mean-based classification

o dNME-phenoName1-VS-phenoName2.bed
       differences in normalized methylation entropies

o DEU-phenoName1-VS-phenoName2.bed
       differential entropy-based classification

o JSD-phenoName1-VS-phenoName2.bed
       Jensen-Shannon distances 

o dESI-phenoName1-VS-phenoName2.bed (if ESIflag = 1) 
       differences in entropic sensitivity indices

o dCAP-phenoName1-VS-phenoName2.bed (if MCflag = 1) 
       differences in channel capacities

o dRDE-phenoName1-VS-phenoName2.bed (if MCflag = 1) 
       differences in relative dissipated energies
