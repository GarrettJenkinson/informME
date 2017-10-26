-----------------------------------------------------------------

# informME

-----------------------------------------------------------------

An information-theoretic pipeline for methylation analysis of WGBS data.

LATEST RELEASE: v0.3.0

-----------------------------------------------------------------

This directory contains the MATLAB, C++ and R software that implements informME, as well as the BASH wrappers and example submission scripts used to run informME on a Sun Grid Engine (SGE) cluster or a SLURM computer cluster. 

informME is an information-theoretic approach to methylation analysis developed in [1] and [2], using the Ising model of statistical physics. This method is applied on BAM files with reads aligned to a reference genome, which are generated during  a whole-genome bisulphite sequencing (WGBS) experiment, and produces genome-wide information about the statistical  properties of methylation in a single-sample and in a differential analysis framework. The resulting information  is stored in multiple MATLAB MAT files for subsequent  processing and is summarized by bedGraph genomic tracks that  can be visualized using a genome browser (such as the UCSC genome browser, see https://genome.ucsc.edu).  

The current implementation of informME has been tested within the following environments:

o Red Hat Enterprise Linux Server Release 6.5 (Santiago) and CentOS release 6.7 (Final)

o Sun Grid Engine OGS/GE 2011.11 and Slurm 17.02.3

o MATLAB R2013b 64-bit and MATLAB R2016b 64-bit  with Bioinformatics and Symbolic Math Toolboxes 

and by using the following tools:

o Bismark Bisulphite Mapper - v0.13.1 (to make BAM files)

o SAMtools - v0.1.19 (to parse BAM files)

o gcc/g++ compiler - v4.4.7-4 and v4.9.2 (to compile MATLAB MEX code)


A. DIRECTORY STRUCTURE
----------------------

informME includes the following directories:

informME/bin - contains softlinks to all executable bash scripts

informME/cluster - contains templates for running informME on SGE and SLURM clusters

informME/src - contains all bash wrappers that include (i) help file, and (ii) comprehensive toy example. Also contains the MATLAB, C++, and R informME source codes

informME/third_party - contains third party MATLAB software used by informME


B. DEPENDENCIES
---------------

It is required to properly configure a MATLAB compiler. See:

https://www.mathworks.com/products/compiler.html 

https://www.mathworks.com/help/matlab/ref/mex.html

for details.

GMP, MPFR, MPREAL, and Eigen are all C++ dependencies that the installation script will install for you as needed.

A working SAMtools installation needs to be on the system path.


C. INSTALLING InformME
----------------------

Run install.sh in the informME directory. During the interactive installing process the user will be asked different questions regarding the locations of default directories. Dependencies, such as GMP, MPFR, MPREAL, and Eigen, will be automatically installed during this process if these are not already available, and informME's C++ MEX code will be compiled.

Note: the following environment variables will be defined in a configuration file stored in ~/.informME/informME.config and then will be accessed throughout multiple points in the informME pipeline:

o REFGENEDIR: directory where the CpGlocationChrX.mat files are stored

o BAMDIR: directory where the BAM files are stored

o INTERDIR: directory where all the intermediate files are stored

o FINALDIR: directory where the output BED files are stored
 
o MATLICENSE: path to MATLAB license

These environment variables can be overwritten through optional arguments when running any part of informME. If these variables are not defined, then the user must pass the corresponding paths as optional arguments.


D. RUNNING informME
-------------------

Run the informME software using the following steps in the indicated order (type the name of the command with no arguments to see help file and instructions in each case).

D.1. REFERENCE GENOME ANALYSIS:
-------------------------------
	
        fastaToCpg.sh [OPTIONS] -- FASTA_FILE

This step analyzes the reference genome and produces a MATLAB MAT file CpGlocationChr#.mat for each chromosome. Each MAT file contains the following information: 

o location of CpG sites 

o CpG density for each CpG site 

o distance between neighboring CpG sites

o location of the last CpG site in the chormosome

o length of chromosome (in base pairs)

NOTE: This step must be completed before proceeding with the next step, but it is only done once per reference genome. 

D.2. METHYLATION DATA MATRIX GENERATION: 
----------------------------------------
	
        getMatrices.sh [OPTIONS]  -- BAM_FILE CHR_NUM TOTAL_PROC

Wrapper that takes a BAM file as input and generates matrices for informME. Each chromosome subdirectory includes the MATLAB file phenoName_matrices.mat. Each file contains the following information for each genomic region, which will be subsequently used in model estimation:

o data matrix with -1,0,1 values for methylation status

o CpG locations broken down by region.

NOTE: This step must be completed before proceeding with the next step. 

D.3. MODEL ESTIMATION & ANALYSIS:
---------------------------------

        informME_run.sh [OPTIONS] -- MAT_FILES PREFIX CHR_NUM TOTAL_PROC

This step estimates the parameters of the Ising probability distribution used to model methylation within equally sized (in base pairs) non-overlapping regions of the genome. Each phenoName_fit.mat file contains the following information for each genomic region used in model estimation: 

o CpG distances

o CpG densities

o estimated alpha, beta, and gamma parameters of the Ising 
model

o initial and transition probabilities of the inhomogeneous 
Markov chain representation of the Ising model

o marginal probabilities at each CpG site

o the log partition function of the estimated Ising model

This is followed by methylation analysis of a given phenotype by computing a number of statistical summaries of the methylation state, including probability distributions of methylation levels, mean methylation levels, and normalized methylation entropies, as well as mean and entropy based classifications. If desired, this step also computes entropic sensitivity indices, as well information-theoretic quantities associated with methylation channels, such as turnover ratios, channel capacities, and relative dissipated energies. Each chromosome subdirectory contains the MATLAB file phenoName_analysis.mat for each phenotypic methylation sample (lungnormal-1, lungcancer-1,  etc.), with each phenoName_analysis.mat file containing the following information for each genomic region used in model estimation: 

o The locations of the CpG sites within the genomic region

o Numbers of CpG sites within the analysis subregions 

o Which analysis subregions are modeled and which are not

o Estimated parameters of Ising model in genomic region

o Methylation level probabilities in modeled subregions

o Coarse methylation level probabilities

o Mean methylation levels

o Normalized methylation entropies

o Entropic sensitivity indices (if ESIflag = 1)

o Turnover ratios (if MCflag = 1)

o Channel capacities (if MCflag = 1)

o Relative dissipated energies (if MCflag = 1)

NOTE: This step must be completed before proceeding with the next step. 

D.4. GENERATE BED FILES FOR SINGLE ANALYSIS:
------------------------------------------

        singleMethAnalysisToBed.sh [OPTIONS] -- PREFIX

This function makes BED files for the methylation analysis results obtained after running informME_run.sh for a given phenotype.

D.5. GENERATE BED FILES FOR DIFFERENTIAL ANALYSIS:
------------------------------------------------

	makeBedsForDiffMethAnalysis.sh [OPTIONS] -- PREFIX_1 PREFIX_2

This function makes BED files for the methylation analysis results obtained after running informME_run.sh for two given phenotypes.

D.6. POST-PROCESSING
------------------

D.6.1. BED TO BW CONVERSION
---------------------------

The user can employ a provided utility to convert BED files generated by informME to much smaller BigWig (BW) files. Run the command:

    ./bed2bw.sh "path" "asy"

where 

- "path" is the path containing the BED files

- "asy" is the assembly of the reference genome used to generate the BED files (for example, "asy" must be set to "hg19" when the Human assembly hg19 is used)

NOTE: For this utility, the following tools must be installed in $PATH: bedtools, bedClip, bedGraphToBigWig, and fetchChromSizes

D.6.2. DMR DETECTION
--------------------

The user can use a provided utility to perform DMR detection using the Jensen-Shannon distance (JSD) based on the method described in [2]. This utility must be run within an R session.

usage (when replicate reference data is available): 

    setwd("/path/to/informME/src/R_src/")
    source("jsDMR.R") 
    runReplicateDMR(refVrefFiles,testVrefFiles,
                    inFolder,outFolder)

where 

- refVrefFiles is a vector of BED file names that contain the JSD values of all pairwise reference comparisons 

- testVrefFiles is a vector of BED file names that contain the JSD values of test/reference comparisons

- inFolder is the directory that contains all JSD files

- outFolder is the directory used to write the results

usage (when no replicate reference data is available) 

    setwd("/path/to/informME/src/R_src/")
    source("jsDMR.R") 
    runNoReplicateDMR(JSDfile,inFolder,outFolder)

where

- JSDfile is the name of a BED file that contains the JSD values of a test/reference comparison

- inFolder is the directory that contains the JSD file

- outFolder is the directory used to write the result

NOTE: For this utility, the following tools must be installed in R: rtracklayer, logitnorm, mixtools.

D.6.3. GENE RANKING
-------------------

The user can use a provided utility to rank all Human genes in the Bioconductor library TxDb.Hsapiens.UCSC.hg19.knownGene using the Jensen-Shannon distance (JSD) based on the method described in [2]. This utility must be run within an R session.

usage (when replicate reference data is available):

    setwd("path/to/informME/src/R_src/")
    source("jsGrank.R")
    rankGenes(refVrefFiles,testVrefFiles,inFolder,outFolder,
              tName,rName)

where 

- refVrefFiles is a vector of BED files that contain the JSD values of a test/reference comparison

- testVrefFiles is a vector of BED files that contain the JSD values of available test/reference comparisons

- inFolder is the directory that contains the JSD files

- outFolder is the directory used to write the result in an .xlsx file

- tName is a string providing a name for the test phenotype

- rName is a string providing a name for the reference phenotype

usage (when no replicate reference data is available):  

    setwd("path/to/informME/src/R_src/")
    source("jsGrank.R")
    rankGenes(c(),testVrefFiles,inFolder,outFolder,
              tName,rName)

where 

- testVrefFiles is a vector of BED files that contain the JSD values of available test/reference comparisons

- inFolder is the directory that contains the JSD files
outFolder is the directory used to write the result in an .xlsx file

- tName is a string providing a name for the test phenotype

- rName is a string providing a name for the reference phenotype

NOTE: For this utility, the following tools must be installed in R: GenomicFeatures, GenomicRanges, Homo.sapiens, rtracklayer, TxDb.Hsapiens.UCSC.hg19.knownGene, XLConnect.


E. OUTPUT FILES
---------------

informME generates the following methylation tracks in the form of BED files.

E.1 SINGLE SAMPLE ANALYSIS
--------------------------

  o MML-PhenoName.bed: mean methylation levels

  o NME-PhenoName.bed: normalized methylation entropy

  o METH-PhenoName.bed: methylation-based classification (non-variable)

  o VAR-PhenoName.bed: methylation-based classification (variable)

  o ENTR-PhenoName.bed: entropy-based classification

  o ESI-PhenoName.bed: entropic sensitivity indices (if ESIflag = 1)

  o TURN-PhenoName.bed: turnover ratios (if MCflag = 1) 

  o CAP-PhenoName.bed: channel capacities(if MCflag = 1)  

  o RDE-PhenoName.bed: relative dissipated energies (if MCflag = 1)

E.2 DIFFERENTIAL ANALYSIS
-------------------------

  o dMML-tPhenoName-VS-rPhenoName.bed: differences in mean methylation levels

  o dNME-tPhenoName-VS-rPhenoName.bed: differences in normalized methylation entropies

  o DMU-tPhenoName-VS-rPhenoName.bed: differential mean-based classification

  o DEU-tPhenoName-VS-rPhenoName.bed: differential entropy-based classification

  o JSD-tPhenoName-VS-rPhenoName.bed: Jensen-Shannon distances

  o dESI-tPhenoName-VS-rPhenoName.bed: differences in entropic sensitivity indices (if ESIflag = 1) 

  o dCAP-tPhenoName-VS-rPhenoName.bed: differences in channel capacities (if MCflag = 1) 

  o dRDE-tPhenoName-VS-rPhenoName.bed: differences in relative dissipated energies (if MCflag = 1) 

E.3 POST PROCESSING
-------------------

In addition, the following files are generated when using the provided post-processing utilities: 

  o DMR-JSD-tPhenoName-VS-rPhenoName.bed: differentially methylated regions (when the jsDMR.R utility is used)

  o gRank-JSD-tName-VS-rName.xlsx: JSD based gene ranking when no replicate reference data is available (when the jsGRank.R utility is used)

  o gRankRRD-JSD-tName-VS-rName.xlsx: JSD based gene ranking when replicate reference data is available (when the jsGRank.R utility is used)


REFERENCES
----------

[1] Jenkinson, G., Pujadas, E., Goutsias, J., and Feinberg, A.P. (2017), Potential energy landscapes indentify the information-theoretic nature of the epigenome, Nature Genetics, 49: 719-729.

[2] Jenkinson, G., Abante, J., Feinberg, A.P., and Goutsias, J. (2017), An information-theoretic approach to the modeling and analysis of whole-genome bisulfite sequencing data, Submitted.


VERSION HISTORY
---------------

v0.3.0 - Code reorganized into a source and bin directory structure. Wrapper shell scripts added to all user callable functions. I/O directories no longer inside code directories. General UI improvements for easier usage.

v0.2.1 - Added R utilities for DMR detection and gene ranking using the Jensen-Shannon distance (JSD). Various documentation improvements.

v0.2.0 - Code reorganized into more specialized directories, streamlined, and general SGE submission scripts provided as a guide. Updated README's. Variable names changed to reflect published notation.

v0.1.0 - Initial release. Code widely tested internally. Code used to create results for ref [1].


LICENCING
---------

All code authored by Garrett Jenkinson or Jordi Abante in informME is licensed under a GPLv3 license; exceptions to GPL licensing are the files contained in the following directories:

informme/third_party/global_optim

informme/third_party/maxent

These files have their own licensing information in their headers. Thanks to Arnold Neumaier and Ali Mohammad-Djafari for their permissions to modify and distribute their software with informME.

-----------------------------------------------------------------

