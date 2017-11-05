-----------------------------------------------------------------

# informME

## An information-theoretic pipeline for methylation analysis of WGBS data.

### LATEST RELEASE: v0.3.0
-----------------------------------------------------------------

This directory contains the MATLAB, C++ and R software that implements informME, as well as the BASH wrappers and example submission scripts used to run informME on a Sun Grid Engine (SGE) cluster or a SLURM computer cluster. 

informME is an information-theoretic approach to methylation analysis developed in [1] and [2], using the Ising model of statistical physics. This method is applied on BAM files with reads aligned to a reference genome, which are generated during  a whole-genome bisulphite sequencing (WGBS) experiment, and produces genome-wide information about the statistical  properties of methylation in a single-sample and in a differential analysis framework. The resulting information  is stored in multiple MATLAB MAT files for subsequent  processing and is summarized by bedGraph genomic tracks that  can be visualized using a genome browser (such as the UCSC genome browser, see https://genome.ucsc.edu).  

The current implementation of informME has been tested within the following environments:

* Red Hat Enterprise Linux Server Release 6.5 (Santiago) and CentOS release 6.7 (Final)

* Sun Grid Engine OGS/GE 2011.11 and Slurm 17.02.3

* MATLAB R2013b 64-bit, MATLAB R2016b 64-bit, and MATLAB R2017a 64-bit with Bioinformatics and Symbolic Math Toolboxes 

and by using the following tools:

* Bismark Bisulphite Mapper - v0.13.1 (to make BAM files)

* SAMtools - v0.1.19 (to parse BAM files)

* gcc/g++ compiler - v4.4.7-4 and v4.9.2 (to compile MATLAB MEX code)


A. DIRECTORY STRUCTURE
----------------------

informME includes the following directories:

1. informME/bin - contains softlinks to all executable bash scripts

2. informME/src - contains all bash wrappers as well as the MATLAB, C++, and R code used by informME

3. informME/third\_party - contains third party MATLAB software used by informME

4. informME/cluster - contains templates for running informME on SGE and SLURM clusters

B. DEPENDENCIES
---------------

It is required to properly configure a MATLAB compiler. See:

* https://www.mathworks.com/products/compiler.html 

* https://www.mathworks.com/help/matlab/ref/mex.html

for details. A working SAMtools installation needs to be on the system path as well. The remaining dependencies (GMP, MPFR, MPREAL, and Eigen) are all C++ dependencies that the installation script will install as needed.

C. INSTALLING InformME
----------------------

Run install.sh in the informME directory. During the interactive installing process the user will be asked different questions regarding the locations of default directories. Dependencies, such as GMP, MPFR, MPREAL, and Eigen, will be automatically installed during this process if these are not already available, and the C++ MEX code of informME will be compiled.

Note 1: the following environment variables will be defined in a configuration file stored in ~/.informME/informME.config and then will be accessed throughout multiple points in the informME pipeline:

* REFGENEDIR: directory where the CpGlocationChrX.mat files will be stored

* BAMDIR: directory where the BAM files are stored

* INTERDIR: directory where all the intermediate files will be stored

* FINALDIR: directory where the output BED files will be stored
 
* MATLICENSE: path to MATLAB license

Note 2: These environment variables can be overwritten through optional arguments when running any part of informME. If these variables are not defined, then the user must pass the corresponding paths as optional arguments.


D. RUNNING informME
-------------------

Run the informME software using the following steps in the indicated order (type the name of the command with no arguments to see the help file with instructions in each case).

D.1. REFERENCE GENOME ANALYSIS:
-------------------------------
	
        fastaToCpg.sh [OPTIONS] -- FASTA_FILE

This step analyzes the reference genome FASTA\_FILE (in FASTA format) and produces a MATLAB MAT file CpGlocationChr#.mat for each chromosome. Each MAT file produced will be stored by default in REFGENEDIR, and it will contain the following information:

* location of CpG sites 

* CpG density for each CpG site 

* distance between neighboring CpG sites

* location of the last CpG site in the chromosome

* length of chromosome (in base pairs)


D.2. METHYLATION DATA MATRIX GENERATION: 
----------------------------------------
	
        getMatrices.sh [OPTIONS]  -- BAM_FILE CHR_NUM

This step takes the BAM file BAM\_FILE as input and generates the methylation data matrix for chromosome number CHR\_NUM. The file BAM\_FILE is expected to be in BAMDIR by default. The output produced by this step is stored by default in a subdirectory in INTERDIR named after the chromosome number CHR\_NUM. The file will preserve the prefix from the file BAM\_FILE and the suffix '\_matrices.mat' will be appended (e.g. if BAM\_FILE is normal\_sample.bam and CHR\_NUM is 10, then the output file will be stored as INTERDIR/chr10/normal\_sample\_matrices.mat by default). The output file profuced will contain, for each genomic region, the following information which will be subsequently used in model estimation:

* data matrix with -1,0,1 values for methylation status

* CpG locations broken down by region.

NOTE: we recommend taking advantage of the array feature available in SGE and SLURM based clusters to submit an individual job for each chromosome.


D.3. MODEL ESTIMATION & ANALYSIS:
---------------------------------

        informME_run.sh [OPTIONS] -- MAT_FILES PHENO CHR_NUM

This step is comprised of two phases. In the first phase, informME learns the parameters of the Ising probability distribution by combining the methylation data matrices provided through the argument MAT\_FILES (comma-separated list) for chromosome number CHR\_NUM. The MAT\_FILES are all expected to be in a subdirectory named after CHR\_NUM in INTERDIR by default. The output generated by this phase will be stored by default in a subdirectory in INTERDIR named after chromosome number CHR\_NUM as well, and will have as prefix PHENO and the suffix '\_fit.mat' will be appended (e.g. if sample\_normal-1,sample\_normal-2,sample\_normal-3 is the list passed as MAT\_FILES, 'normal' is the PHENO passed, and CHR\_NUM is 10, then the output will be stored as INTERDIR/chr10/normal\_fit.mat). The output file will contain the following information:

* CpG distances

* CpG densities

* estimated alpha, beta, and gamma parameters of the Ising model

* initial and transition probabilities of the inhomogeneous Markov chain representation of the Ising model

* marginal probabilities at each CpG site

* the log partition function of the estimated Ising model

The second phase of this step consists in analyzing the model learned by computing a number of statistical summaries of the methylation state, including probability distributions of methylation levels, mean methylation levels, and normalized methylation entropies, as well as mean and entropy based classifications. If desired, this step also computes entropic sensitivity indices, as well information-theoretic quantities associated with methylation channels, such as turnover ratios, channel capacities, and relative dissipated energies. The output of this second phase will be stored in the same directory as the otuput of the first phase, and will have the same prefix as well. However, the suffix in this case will be '\_analysis.mat' (e.g. following the previous example, the path of the output file of this second phase will be INTERDIR/chr10/normal\_analysis.mat). The output file will contain the following information:

* The locations of the CpG sites within the genomic region

* Numbers of CpG sites within the analysis subregions 

* Which analysis subregions are modeled and which are not

* Estimated parameters of Ising model in genomic region

* Methylation level probabilities in modeled subregions

* Coarse methylation level probabilities

* Mean methylation levels

* Normalized methylation entropies

* Entropic sensitivity indices (if ESIflag = 1)

* Turnover ratios (if MCflag = 1)

* Channel capacities (if MCflag = 1)

* Relative dissipated energies (if MCflag = 1)

NOTE: we recommend taking advantage of the array feature available in SGE and SLURM based clusters to submit an individual job for each chromosome.

D.4. GENERATE BED FILES FOR SINGLE ANALYSIS:
------------------------------------------

        singleMethAnalysisToBed.sh [OPTIONS] -- PHENO

This function makes BED files for the methylation analysis results obtained after running informME\_run.sh for a given phenotype PHENO. The input files (analysis file) are expected to have the path INTERDIR/chr#/PHENO\_analysis.mat by default. The output files will be stored by default in FINALDIR and will have the following names and content:

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

	makeBedsForDiffMethAnalysis.sh [OPTIONS] -- PHENO1 PHENO2

This function makes BED files for the methylation analysis results obtained after running informME\_run.sh for two given phenotypes PHENO1 and PHENO2. The input files (both analysis files) are expected to have the path INTERDIR/chr#/PHENO1\_analysis.mat and INTERDIR/chr#/PHENO2\_analysis.mat by default. The output files will be stored by default in FINALDIR and will have the following names and content:

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

The user can employ a provided utility to convert BED files generated by informME to much smaller BigWig (BW) files. Run the command:

    ./bed2bw.sh "path" "asy"

where 

* "path" is the path containing the BED files

* "asy" is the assembly of the reference genome used to generate the BED files (for example, "asy" must be set to "hg19" when the Human assembly hg19 is used)

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

* refVrefFiles is a vector of BED file names that contain the JSD values of all pairwise reference comparisons 

* testVrefFiles is a vector of BED file names that contain the JSD values of test/reference comparisons

* inFolder is the directory that contains all JSD files

* outFolder is the directory used to write the results

usage (when no replicate reference data is available) 

    setwd("/path/to/informME/src/R_src/")
    source("jsDMR.R") 
    runNoReplicateDMR(JSDfile,inFolder,outFolder)

where

* JSDfile is the name of a BED file that contains the JSD values of a test/reference comparison

* inFolder is the directory that contains the JSD file

* outFolder is the directory used to write the result

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

* refVrefFiles is a vector of BED files that contain the JSD values of a test/reference comparison

* testVrefFiles is a vector of BED files that contain the JSD values of available test/reference comparisons

* inFolder is the directory that contains the JSD files

* outFolder is the directory used to write the result in an .xlsx file

* tName is a string providing a name for the test phenotype

* rName is a string providing a name for the reference phenotype

usage (when no replicate reference data is available):  

    setwd("path/to/informME/src/R_src/")
    source("jsGrank.R")
    rankGenes(c(),testVrefFiles,inFolder,outFolder,
              tName,rName)

where 

* testVrefFiles is a vector of BED files that contain the JSD values of available test/reference comparisons

* inFolder is the directory that contains the JSD files
outFolder is the directory used to write the result in an .xlsx file

* tName is a string providing a name for the test phenotype

* rName is a string providing a name for the reference phenotype

NOTE: For this utility, the following tools must be installed in R: GenomicFeatures, GenomicRanges, Homo.sapiens, rtracklayer, TxDb.Hsapiens.UCSC.hg19.knownGene, XLConnect.


E. OUTPUT FILES
---------------

informME generates the following methylation tracks in the form of BED files.

E.1 SINGLE SAMPLE ANALYSIS
--------------------------

* MML-PhenoName.bed: mean methylation levels

* NME-PhenoName.bed: normalized methylation entropy

* METH-PhenoName.bed: methylation-based classification (non-variable)

* VAR-PhenoName.bed: methylation-based classification (variable)

* ENTR-PhenoName.bed: entropy-based classification

* ESI-PhenoName.bed: entropic sensitivity indices (if ESIflag = 1)

* TURN-PhenoName.bed: turnover ratios (if MCflag = 1) 

* CAP-PhenoName.bed: channel capacities(if MCflag = 1)  

* RDE-PhenoName.bed: relative dissipated energies (if MCflag = 1)

E.2 DIFFERENTIAL ANALYSIS
-------------------------

* dMML-tPhenoName-VS-rPhenoName.bed: differences in mean methylation levels

* dNME-tPhenoName-VS-rPhenoName.bed: differences in normalized methylation entropies

* DMU-tPhenoName-VS-rPhenoName.bed: differential mean-based classification

* DEU-tPhenoName-VS-rPhenoName.bed: differential entropy-based classification

* JSD-tPhenoName-VS-rPhenoName.bed: Jensen-Shannon distances

* dESI-tPhenoName-VS-rPhenoName.bed: differences in entropic sensitivity indices (if ESIflag = 1) 

* dCAP-tPhenoName-VS-rPhenoName.bed: differences in channel capacities (if MCflag = 1) 

* dRDE-tPhenoName-VS-rPhenoName.bed: differences in relative dissipated energies (if MCflag = 1) 

E.3 POST PROCESSING
-------------------

In addition, the following files are generated when using the provided post-processing utilities: 

* DMR-JSD-tPhenoName-VS-rPhenoName.bed: differentially methylated regions (when the jsDMR.R utility is used)

* gRank-JSD-tName-VS-rName.xlsx: JSD based gene ranking when no replicate reference data is available (when the jsGRank.R utility is used)

* gRankRRD-JSD-tName-VS-rName.xlsx: JSD based gene ranking when replicate reference data is available (when the jsGRank.R utility is used)


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

* informme/third\_party/global\_optim

* informme/third\_party/maxent

These files have their own licensing information in their headers. Thanks to Arnold Neumaier and Ali Mohammad-Djafari for their permissions to modify and distribute their software with informME.

-----------------------------------------------------------------

