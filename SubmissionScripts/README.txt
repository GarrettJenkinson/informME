Last Modified: 04/22/17

----------------------------------------------------------------
SUBDIRECTORY SubmissionScripts
----------------------------------------------------------------

This directory contains the submission scripts used to run 
informME on a computer cluster using the Sun Grid Engine 
(SGE) environment. These files can be used as templates to 
develop scripts that can allow the user to run informME 
on other grid engine parallel enviroments. This directory 
is organized using the following subdirectories:

MexSource
---------

The following script should be run within this subdirectory 
in order to compile and install informME:

./runCompileMEX.sh

** run only once to compile/install the C++ MEX code 
   required by informME

ParseBAMfile
------------

The following script must be run first within this 
subdirectory in order to analyze the reference genome 
(for example, for the reference genome FASTA file 
"Homo_sapiens_assembly19.fasta" and for the species 
"Human"):

./runGenomeAnalysis.sh "Homo_sapiens_assembly19.fasta" "Human" 

** run only once per reference genome

Subsequently, the following script must be run within 
this subdirectory for each BAM file (for example, for 
using file lungnormal-1.bam for species "Human" and 
for trimming 15 bases from the first read in a read 
pair and 17 bases from the second read in the pair):

./runDataMatrixGeneration.sh "lungnormal-1" "Human" "[15,17]"

** run for each BAM file to generate methylation data 
   matrices 

Modeling
--------

The following script must be run within this subdirectory 
to build statistical models of data from BAM files for 
which the "ParseBAMfile" code has already been run (for 
example, for modeling two parsed BAM files, 
lungnormal-1.bam and lungnormal-2.bam that provide 
methylation data for a normal lung phenotype in the 
species "Human"):

./runModelEstimation.sh "{'lungnormal-1','lungnormal-2'}" "lungnormal" "Human"

** run for each model to be built from a set of parsed BAM 
   files

SingleAnalysis
--------------

The following script must be run within this subdirectory 
to analyze a single model that was built using the 
"Modeling" code (for example, for a model built from a 
normal lung phenotype in the species "Human"):

./runSingleMethAnalysis.sh "lungnormal-1" "Human"

** run for each model to generate single (inter-sample) 
   methylation analysis results

DiffAnalysis
------------

The following script must be run within this subdirectory
to produce a differential analysis of two models that 
have been analyzed using the "SingleAnalysis" code (for 
example, for a test and reference model, which were 
respectively built from a lung cancer and lung normal 
phenotype in the species "Human"):

./runDiffMethAnalysis.sh "lungcancer-1" "lungnormal-1" "Human"

** run for each pair of a test and a reference model to 
   generate differential methylation results

STDouts
-------

Contains standard output text files for troubleshooting.
