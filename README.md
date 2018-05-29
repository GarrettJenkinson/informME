-----------------------------------------------------------------

# informME

## An information-theoretic pipeline for methylation analysis of WGBS data

### LATEST RELEASE: v0.3.1
-----------------------------------------------------------------

This directory contains the MATLAB, C++ and R software that implements informME, as well as the BASH wrappers and example submission scripts used to run informME on a Sun Grid Engine (SGE) cluster or a SLURM computer cluster. 

informME is an information-theoretic approach to methylation analysis developed in [1] and [2], using the Ising model of statistical physics. This method is applied on BAM files with reads aligned to a reference genome, which are generated during  a whole-genome bisulphite sequencing (WGBS) experiment, and produces genome-wide information about the statistical  properties of methylation in a single-sample and in a differential analysis framework. The resulting information  is stored in multiple MATLAB MAT files for subsequent  processing and is summarized by bedGraph genomic tracks that  can be visualized using a genome browser (such as the UCSC genome browser, see https://genome.ucsc.edu).  

The current implementation of informME has been tested within the following environments:

* Red Hat Enterprise Linux Server Release 6.5 (Santiago), CentOS release 6.7 (Final) and Ubuntu 16.04 (xenial)

* Sun Grid Engine OGS/GE 2011.11 and Slurm 17.02.3

* MATLAB R2013b 64-bit, MATLAB R2016b 64-bit, MATLAB R2017a 64-bit, and MATLAB R2018a 64-bit with Bioinformatics and Symbolic Math Toolboxes 

and by using the following tools:

* Bismark Bisulphite Mapper - v0.13.1 (to make BAM files)

* SAMtools - v0.1.19 (to parse BAM files)

* gcc/g++ compiler - v4.4.7-4, v4.9.2, v6.3.0 (to compile MATLAB MEX code, version must be compatible with corresponding MATLAB release, see section B below)


A. DIRECTORY STRUCTURE
----------------------

informME includes the following directories:

1. informME/bin - contains softlinks to all executable bash scripts

2. informME/src - contains all bash wrappers as well as the MATLAB, C++, and R code used by informME

3. informME/third\_party - contains third party MATLAB software used by informME

4. informME/cluster - contains templates for running informME on SGE and SLURM clusters

B. DEPENDENCIES
---------------

It is required to properly configure a MATLAB C++ mex compiler. For details, see:

* https://www.mathworks.com/help/matlab/ref/mex.html

A working SAMtools installation needs to be on the system path as well. The remaining dependencies (GMP, MPFR, MPREAL, and Eigen) are all C++ dependencies that the installation script will install as needed.

NOTE: The required "mex compiler" is different from the "Matlab compiler SDK". InformME requires the former and not the latter.

C. INSTALLING InformME
----------------------

Run install.sh in the informME directory. During the interactive installation process the user will be asked different questions regarding the locations of default directories. Dependencies, such as GMP, MPFR, MPREAL, and Eigen, will be automatically installed during this process if these are not already available, and the C++ MEX code of informME will be compiled.

Note 1: the following environment variables will optionally be defined in a configuration file stored in ~/.informME/informME.config and then will be accessed throughout multiple points in the informME pipeline:

* REFGENEDIR: directory where the CpGlocationChr#.mat files will be stored

* BAMDIR: directory where the BAM files are stored

* SCRATCHDIR: directory where all temporary files will be stored

* INTERDIR: directory where all intermediate files will be stored

* FINALDIR: directory where the output BED files will be stored
 
* MATLICENSE: path to MATLAB license

Note 2: The previous environment variables can be overwritten through optional arguments when running any part of informME. If these variables are not defined, then the user must pass the corresponding paths as optional arguments.


D. RUNNING informME
-------------------

Run the informME software using the following steps in the indicated order (type the name of the command with no arguments to see the help file with instructions in each case).

D.1. REFERENCE GENOME ANALYSIS:
-------------------------------
	
        fastaToCpg.sh [OPTIONS] -- FASTA_FILE

This step analyzes the reference genome FASTA\_FILE (in FASTA format) and produces a MATLAB MAT file CpGlocationChr#.mat for each chromosome, which is stored by default in REFGENEDIR, and contains the following information:

* location of CpG sites 

* CpG density for each CpG site 

* distance between neighboring CpG sites

* location of the last CpG site in the chromosome

* length of chromosome (in base pairs)

NOTE1: This step only needs to be completed one time for a given reference genome. Start analyzing samples at step D.2 if you have previously completed step D.1 for your sample's reference genome.

NOTE2: At this time the statistical model of informME has been designed to work only with autosomes, and so the informME software will not model mitochondrial chromosomes, lambda spike-ins, partial contigs, sex chromosomes, et cetera. Also the reference fasta file to which bam files have been aligned is assumed to be sorted so that the somatic chromosomes come first and in the usual order: chr1,chr2,...,chrN. 

D.2. METHYLATION DATA MATRIX GENERATION: 
----------------------------------------
	
        getMatrices.sh [OPTIONS]  -- BAM_FILE CHR_NUM

This step takes the BAM file BAM\_FILE as input and generates the methylation data matrix for chromosome number CHR\_NUM. By default, the file BAM\_FILE and its associated index file (with extension .bai) is expected to be in BAMDIR, and the output file produced by this step is stored in a subdirectory in INTERDIR named after the chromosome number CHR\_NUM. The output file preserves the prefix from the file BAM\_FILE and the suffix '\_matrices.mat' is appended to it (e.g. if BAM\_FILE is normal\_sample.bam and CHR\_NUM is 10, then the output file is saved as INTERDIR/chr10/normal\_sample\_matrices.mat). The file produced contains the following information for each genomic region, which is subsequently used for model estimation:

* data matrix with -1,0,1 values for methylation status

* CpG locations broken down by region

NOTE1: We recommend taking advantage of the array feature available in SGE and SLURM based clusters to submit an individual job for each chromosome.

NOTE2: See reference [1], "Online Methods: Quality control and alignment" for our suggested preprocessing steps when generating a sorted, indexed, deduplicated BAM file to input to informME.  

D.3. MODEL ESTIMATION & ANALYSIS:
---------------------------------

        informME_run.sh [OPTIONS] -- MAT_FILES PHENO CHR_NUM

This step is comprised of two phases. During the first phase, informME learns the parameters of the Ising probability distribution by combining the methylation data matrices provided through the argument MAT\_FILES (comma-separated list) for chromosome number CHR\_NUM. By default, the MAT\_FILES are expected to be in a subdirectory named after CHR\_NUM in INTERDIR. The output generated during this phase is also stored in a subdirectory in INTERDIR named after chromosome number CHR\_NUM. The output file has as prefix PHENO and the suffix '\_fit.mat' appended to it (e.g. if 'normal' is the PHENO, and CHR\_NUM is 10, then the output is stored as INTERDIR/chr10/normal\_fit.mat). The file produced contains the following information:

* CpG distances

* CpG densities

* estimated alpha, beta, and gamma parameters of the Ising model

* initial and transition probabilities of the inhomogeneous Markov chain representation of the Ising model

* marginal probabilities at each CpG site

* the log partition function of the estimated Ising model

The second phase of this step consists in analyzing the model learned by computing a number of statistical summaries of the methylation state, including probability distributions of methylation levels, mean methylation levels, and normalized methylation entropies, as well as mean and entropy based classifications. If desired, this step also computes entropic sensitivity indices, as well information-theoretic quantities associated with methylation channels, such as turnover ratios, channel capacities, and relative dissipated energies. The output generated during this phase is stored in the same directory as the output generated during the first phase, using the same prefix as before. However, the suffix is now '\_analysis.mat' (e.g. following the previous example, the output file of this phase is stored as INTERDIR/chr10/normal\_analysis.mat). This file contains the following information:

* the locations of the CpG sites within the genomic region

* numbers of CpG sites within the analysis subregions 

* which analysis subregions are modeled and which are not

* estimated parameters of Ising model in genomic region

* methylation level probabilities in modeled subregions

* coarse methylation level probabilities

* mean methylation levels

* normalized methylation entropies

* entropic sensitivity indices (if ESIflag = 1)

* turnover ratios (if MCflag = 1)

* channel capacities (if MCflag = 1)

* relative dissipated energies (if MCflag = 1)

NOTE: We recommend taking advantage of the array feature available in SGE and SLURM based clusters to submit an individual job for each chromosome.

D.4. GENERATE BED FILES FOR SINGLE ANALYSIS:
------------------------------------------

        singleMethAnalysisToBed.sh [OPTIONS] -- PHENO

This function makes BED files from the methylation analysis results obtained after running informME\_run.sh for a given phenotype PHENO. By default, the input file (analysis file) is expected to be located in INTERDIR/chr#/PHENO\_analysis.mat. In addition, the output files are stored in FINALDIR and have the following names and content:

* MML-PHENO.bed: mean methylation levels

* NME-PHENO.bed: normalized methylation entropy
    
* METH-PHENO.bed: methylation-based classification (non-variable)
    
* VAR-PHENO.bed: methylation-based classification (variable)
    
* ENTR-PHENO.bed: entropy-based classification
    
* ESI-PHENO.bed (if ESIflag = 1): entropic sensitivity indices
    
* TURN-PHENO.bed (if MCflag = 1): turnover ratios
    
* CAP-PHENO.bed (if MCflag = 1): channel capacities
    
* RDE-PHENO.bed (if MCflag = 1): relative dissipated energies  


D.5. GENERATE BED FILES FOR DIFFERENTIAL ANALYSIS:
------------------------------------------------

	diffMethAnalysisToBed.sh [OPTIONS] -- PHENO1 PHENO2

This function makes BED files for the differential methylation analysis results obtained after running informME\_run.sh for two given phenotypes PHENO1 and PHENO2. By default, the input files (both analysis files) are expected to be located in INTERDIR/chr#/PHENO1\_analysis.mat and INTERDIR/chr#/PHENO2\_analysis.mat respectively. In addition, the output files are stored in FINALDIR and have the following names and content:

* dMML-PHENO1-VS-PHENO2.bed: differences in mean methylation levels
       
* DMU-PHENO1-VS-PHENO2.bed: differential mean-based classification
       
* dNME-PHENO1-VS-PHENO2.bed: differences in normalized methylation entropies
       
* DEU-PHENO1-VS-PHENO2.bed: differential entropy-based classification
       
* JSD-PHENO1-VS-PHENO2.bed: Jensen-Shannon distances
       
* dESI-PHENO1-VS-PHENO2.bed (if ESIflag = 1): differences in entropic sensitivity indices
       
* dCAP-PHENO1-VS-PHENO2.bed (if MCflag = 1): differences in channel capacities
       
* dRDE-PHENO1-VS-PHENO2.bed (if MCflag = 1): differences in relative dissipated energies


D.6. POST-PROCESSING
------------------

D.6.1. BED TO BW CONVERSION
---------------------------

The user can employ a provided utility to convert BED files generated by informME to much smaller BigWig (BW) files. This can be done by running the command:

    ./bed2bw.sh "path" "asy"

where 

* "path" is the path containing the BED files

* "asy" is the assembly of the reference genome used to generate the BED files (for example, "asy" must be set to "hg19" when the Human assembly hg19 is used)

NOTE 1: Type bed2bw.sh to get more information about this utility

NOTE 2: For this utility, the following tools must be installed in $PATH: bedtools, bedClip, bedGraphToBigWig, and fetchChromSizes

D.6.2. DMR DETECTION
--------------------

The user can use a provided utility to perform DMR detection using the Jensen-Shannon distance (JSD) based on the method described in [2]. This utility must be run within an R session.

usage (when replicate reference data is available): 

    setwd("/path/to/informME/src/R_src/")
    source("jsDMR.R") 
    runReplicateDMR(refVrefFiles,testVrefFiles,inFolder,outFolder)

where 

* refVrefFiles is a vector of BIGWIG file names that contain the JSD values of all pairwise reference comparisons 

* testVrefFiles is a vector of BIGWIG file names that contain the JSD values of test/reference comparisons

* inFolder is the directory that contains all JSD files

* outFolder is the directory used to write the results

usage (when no replicate reference data is available) 

    setwd("/path/to/informME/src/R_src/")
    source("jsDMR.R") 
    runNoReplicateDMR(JSDfile,inFolder,outFolder)

where

* JSDfile is the name of a BIGWIG file that contains the JSD values of a test/reference comparison

* inFolder is the directory that contains the JSD file

* outFolder is the directory used to write the result

NOTE 1: More information about these utilities can be found in informME/src/R\_src/postprocess/README.txt

NOTE 2: For this utilities, the following tools must be installed in R: rtracklayer, logitnorm, mixtools, annotatr, and Homo.sapiens

D.6.3. GENE RANKING
-------------------

The user can use a provided utility to rank all Human genes in the Bioconductor library TxDb.Hsapiens.UCSC.hg19.knownGene using the Jensen-Shannon distance (JSD) based on the method described in [2]. This utility must be run within an R session.

usage (when replicate reference data is available):

    setwd("path/to/informME/src/R_src/")
    source("jsGrank.R")
    rankGenes(refVrefFiles,testVrefFiles,inFolder,outFolder,tName,rName)

where 

* refVrefFiles is a vector of BIGWIG files that contain the JSD values of a test/reference comparison

* testVrefFiles is a vector of BIGWIG files that contain the JSD values of available test/reference comparisons

* inFolder is the directory that contains the JSD files

* outFolder is the directory used to write the result in an .xlsx file

* tName is a string providing a name for the test phenotype

* rName is a string providing a name for the reference phenotype

In this case, the function generates the file gRank-JSD-tName-VS-rName.xlsx.

usage (when no replicate reference data is available):  

    setwd("path/to/informME/src/R_src/")
    source("jsGrank.R")
    rankGenes(c(),testVrefFiles,inFolder,outFolder,tName,rName)

where 

* testVrefFiles is a vector of BIGWIG files that contain the JSD values of available test/reference comparisons

* inFolder is the directory that contains the JSD files
outFolder is the directory used to write the result in an .xlsx file

* tName is a string providing a name for the test phenotype

* rName is a string providing a name for the reference phenotype

In this case, the function generates the file gRankRRD-JSD-tName-VS-rName.xlsx.

NOTE 1: More information about this utility can be found in informME/src/R\_src/postprocess/README.txt

NOTE 2: For this utility, the following tools must be installed in R: GenomicFeatures, GenomicRanges, Homo.sapiens, rtracklayer, TxDb.Hsapiens.UCSC.hg19.knownGene, XLConnect

REFERENCES
----------

[1] Jenkinson, G., Pujadas, E., Goutsias, J., and Feinberg, A.P. (2017), Potential energy landscapes indentify the information-theoretic nature of the epigenome, Nature Genetics, 49: 719-729.

[2] Jenkinson, G., Abante, J., Feinberg, A.P., and Goutsias, J. (2018), An information-theoretic approach to the modeling and analysis of whole-genome bisulfite sequencing data, BMC Bioinformatics, 19:87, https://doi.org/10.1186/s12859-018-2086-5.

VERSION HISTORY
---------------

v0.3.1 - More informative error messages. Fix UI to expose ability to model single-ended sequencing reads. Documentation improvements.

v0.3.0 - Code reorganized into a source and bin directory structure. Wrapper shell scripts added to all user callable functions. I/O directories no longer inside code directories. General UI improvements for easier usage.

v0.2.2 -    Minor updates in gene ranking and small updates getting ready for the next major reorganization release.

v0.2.1 -   	Added R utilities for DMR detection and gene ranking using the Jensen-Shannon distance (JSD). Various documentation improvements.

v0.2.0 - Code reorganized into more specialized directories, streamlined, and general SGE submission scripts provided as a guide. Updated README's. Variable names changed to reflect published notation.

v0.1.0 - Initial release. Code widely tested internally. Code used to create results for ref [1].


LICENCING
---------

All code authored by Garrett Jenkinson or Jordi Abante in informME is licensed under a GPLv3 license; exceptions to GPL licensing are the files contained in the following directories:

* informME/third\_party/global\_optim

* informME/third\_party/maxent

These files have their own licensing information in their headers. Thanks to Arnold Neumaier and Ali Mohammad-Djafari for their permissions to modify and distribute their software with informME.

-----------------------------------------------------------------

