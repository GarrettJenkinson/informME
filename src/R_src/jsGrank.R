#
# informME: An information-theoretic pipeline for WGBS data
# Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# or see <http://www.gnu.org/licenses/>.
#
#########################################################################
#
# This is an R script that ranks all Human genes in the
# Bioconductor library TxDb.Hsapiens.UCSC.hg19.knownGene using
# the Jensen-Shannon distance (JSD) based on the method described
# in [1]. It should be run within an R session.
#
#  usage (replicate reference data is available):
#
#   setwd("path/to/informME/PostProcess/")
#   source("jsGrank.R")
#   rankGenes(refVrefFiles,testVrefFiles,inFolder,outFolder,
#             tName,rName)
#
#   # refVrefFiles is a vector of BED files that contain the
#   # JSD values of a test/reference comparison.
#   # For example: if
#   #
#   # JSD-lungnormal-1-VS-lungnormal-2.bed
#   # JSD-lungcancer-3-VS-lungnormal-1.bed
#   # JSD-lungnormal-3-VS-lungnormal-2.bed
#   #
#   # are available, then set
#   #
#   # textVrefFiles <- c("JSD-lungnormal-1-VS-lungnormal-2.bed",
#   #                    "JSD-lungnormal-3-VS-lungnormal-1.bed",
#   #                    "JSD-lungnormal-3-VS-lungnormal-2.bed")
#   #
#   # testVrefFiles is a vector of BED files that contain the
#   # JSD values of available test/reference comparisons.
#   # For example: if
#   #
#   # JSD-lungcancer-1-VS-lungnormal-1.bed
#   # JSD-lungcancer-2-VS-lungnormal-2.bed
#   # JSD-lungcancer-3-VS-lungnormal-3.bed
#   #
#   # are available, then set
#   #
#   # textVrefFiles <- c("JSD-lungcancer-1-VS-lungnormal-1.bed",
#   #                    "JSD-lungcancer-2-VS-lungnormal-2.bed",
#   #                    "JSD-lungcancer-3-VS-lungnormal-3.bed")
#   #
#   # inFolder is the directory that contains the JSD files
#   # outFolder is the directory used to write the result
#   # (a .xlsx file).
#   #
#   # For example:
#   #
#   # inFolder  <- "/path/to/in-folder/"
#   # outFolder <- "/path/to/out-folder/"
#   #
#   # tName and rName are strings providing names for the
#   # test and reference phenotypes.
#   #
#   # For example:
#   #
#   # tName <- "lungcancer"
#   # rName <- "lungnormal"
#
#  usage (no replicate reference data is available):
#
#   setwd("path/to/informME/PostProcess/")
#   source("jsGrank.R")
#   rankGenes(c(),testVrefFiles,inFolder,outFolder,
#             tName,rName)
#
#   # testVrefFiles is a vector of BED files that contain the
#   # JSD values of available test/reference comparisons.
#   # For example: if
#   #
#   # JSD-lungcancer-1-VS-lungnormal-1.bed
#   # JSD-lungcancer-2-VS-lungnormal-2.bed
#   # JSD-lungcancer-3-VS-lungnormal-3.bed
#   #
#   # are available, then set
#   #
#   # textVrefFiles <- c("JSD-lungcancer-1-VS-lungnormal-1.bed",
#   #                    "JSD-lungcancer-2-VS-lungnormal-2.bed",
#   #                    "JSD-lungcancer-3-VS-lungnormal-3.bed")
#   #
#   # inFolder is the directory that contains the JSD files
#   # outFolder is the directory used to write the result
#   # (a .xlsx file).
#   #
#   # For example:
#   #
#   # inFolder  <- "/path/to/in-folder/"
#   # outFolder <- "/path/to/out-folder/"
#   #
#   # tName and rName are strings providing names for the
#   # test and reference phenotypes.
#   #
#   # For example:
#   #
#   # tName <- "lungcancer"
#   # rName <- "lungnormal"
#
#
# REQUIRED PACKAGES:
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(gamlss))
suppressMessages(library("Homo.sapiens"))

# A bug in package means this needs to be run in top environment
gen.Family("SST", "logit")

# Define function to rank gene proms. When no replicate reference data is available
# pass empty vector c() to refVrefFiles.
rankGenes <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,
                      tName="test",rName="ref",
                      GUsize=150,bpBeforeTSS=2000,bpPastTSS=2000,
                      chrsOfInterest=paste("chr",1:22,sep=""),minGUs=10){

  if(length(refVrefFiles)==0){ # Use JSD in promoter to rank the genes.
    result <- rankGenePromsJSD(testVrefFiles,inFolder,outFolder,tName,rName,
                               GUsize=GUsize,bpBeforeTSS=bpBeforeTSS,bpPastTSS=bpPastTSS,
                               chrsOfInterest=chrsOfInterest,minGUs=minGUs)
    return(result)
  }
  # Else use empirical null distribution to rank the genes based on proms and bodies.
  result <- rankGeneBodiesAndProms(refVrefFiles,testVrefFiles,inFolder,outFolder,tName,rName,
                                   chrsOfInterest = chrsOfInterest,
                                   bpPastTSS=bpPastTSS,bpBeforeTSS=bpPastTSS)
} # end rankGenes function

# Define function that computes the sum of SQS values (stored in numvar)
# within each DMR (stored in bins) - adapted from previously available
# documentation of the tileGenome function in GenomicRanges.
binnedSumms <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
} # end binnedSumms function

# Define function that computes the average of JSD values (stored in numvar)
# within each genomic feature (stored in bins) - adapted from previously available
# documentation of the tileGenome function in GenomicRanges.
binnedAverage <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         viewMeans(views)
                       })
  new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
} # end binnedAverage function

# Define function to rank genes when no replicate reference data is available.
rankGenePromsJSD <- function(testVrefFiles,inFolder,outFolder,tName="test",rName="ref",
                             GUsize=150,bpBeforeTSS=2000,bpPastTSS=2000,
                             chrsOfInterest=paste("chr",1:22,sep=""),minGUs=10,
                             doRankProd=FALSE){

  # Check folders for trailing slash, add if missing
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }

  # Define promoter regions.
  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  genesAll<-genes(txdb)
  genes<-keepSeqlevels(genesAll,chrsOfInterest,pruning.mode="coarse")
  sortSeqlevels(genes)
  proms<-promoters(genes,upstream=bpBeforeTSS,downstream=bpPastTSS)

  # Count number of test/reference comparisons.
  numAltComp  <- length(testVrefFiles)

  # Initialize lists of genomic ranges.
  altGRs  <- list()

  # Load files into genomic range lists and arrange the seqlevels.
  for(ind in 1:numAltComp){
    if(substr(testVrefFiles[ind],nchar(testVrefFiles[ind]),nchar(testVrefFiles[ind]))%in% c("w","g","W","G")){
      altGRs[[ind]] <- import.bw(file.path(inFolder,testVrefFiles[ind],fsep=""))
    }else{
      altGRs[[ind]] <- import.bedGraph(file.path(inFolder,testVrefFiles[ind],fsep=""))
    }
    genome(altGRs[[ind]]) <- "hg19"
    sortSeqlevels(altGRs[[ind]])
    altGRs[[ind]] <- keepSeqlevels(altGRs[[ind]],chrsOfInterest,pruning.mode="coarse")
  }

  # Compute the Fisher sum and the degree of freedom used in meta analysis.
  for(ind in 1:numAltComp){
    # Compute run length encodings of tracks for efficient binnedSumms algorithm.
    score <- coverage(altGRs[[ind]],weight="score")
    cov <- coverage(altGRs[[ind]])

    # Compute sum over promoters.
    proms<-binnedSumms(proms,score,paste("score",ind,sep=""))
    proms<-binnedSumms(proms,cov,paste("cov",ind,sep=""))
    mcols(proms)[[paste("score",ind,sep="")]] <- mcols(proms)[[paste("score",ind,sep="")]]/GUsize
    mcols(proms)[[paste("cov",ind,sep="")]] <- ceiling(mcols(proms)[[paste("cov",ind,sep="")]]/GUsize)

    # Remove insufficient data.
    mcols(proms)[[paste("score",ind,sep="")]][mcols(proms)[[paste("cov",ind,sep="")]]< minGUs] <- NA
  }


  # Generate ranked gene lists.
  rankingTabs <- list()
  for (ind in 1:numAltComp){
    sortedIndices <- order(mcols(proms)[[paste("score",ind,sep="")]],decreasing = TRUE)
    sortedOutput <- AnnotationDbi::select(Homo.sapiens, key=as.character(proms[sortedIndices]$gene_id), keytype="ENTREZID", columns=c("SYMBOL"))
    sortedOutput <- sortedOutput[(!duplicated(sortedOutput$SYMBOL))&(!is.na(sortedOutput$SYMBOL)),,drop=FALSE]
    sortedOutput[[paste("rank",ind,sep="")]] <- 1:length(sortedOutput$SYMBOL)
    rankingTabs[[ind]] <- sortedOutput
    rm(list=c("sortedIndices","sortedOutput"))
  }
  rm(list="proms")

  # Generate merged table.
  df <- rankingTabs[[1]]
  if (numAltComp>1){
    for (ind in 2:numAltComp){
      df <- merge(df,rankingTabs[[ind]])
    }
  }

  # Generate rank product.
  df$rankProd <- 1.0
  for(ind in 1:numAltComp){
    df$rankProd <- as.numeric(df$rankProd)*as.numeric(df[[paste("rank",ind,sep="")]])
  }

  # Sort by rankproduct.
  sortedIndices <- order(df$rankProd,decreasing=FALSE)
  df <- df[sortedIndices,]

  # Compute pvalues if flag is set.
  # Note that this requires the package rankprodbounds.R and
  # an independence assumption that is not always valid.
  if (numAltComp>1 & doRankProd){
    # Get HesEisBre-14 algorithm for merging ranked lists
    source(paste(inFolder,"rankprodbounds.R",sep=""))

    # Compute pValue.
    df$pVals<-rankprodbounds(df$rankprod,length(df$rankprod),numAltComp,Delta = 'geometric')
  }

  # Reformat table to handle one or more comparisons.
  if (numAltComp>1){
    df <- subset(df,select = -TXID)
    df$rankProd <- 1:length(df$rankProd)
  } else{
    df$rank <- df$rank1
    df <- subset(df,select = c(-TXID,-rankProd,-rank1))
  }

  # Write results to output.
  write.table(df,file=paste(outFolder,"gRankProms-JSD-",tName,"-VS-",rName,".tsv",sep=""),
              row.names = FALSE,quote = FALSE,sep="\t")

  # Return value.
  df
}# End rankGenePromsJSD.

#make the empirical p-value function
constructNullPvals <- function(nullValues){
  # make null dist
  Fn <- ecdf(nullValues)

  # return function to compute p-values
  function(x) 1-Fn(x)
} # end constructNullPvals function


# Define Function to remove promoter region from gene bodies
removePromRegion <- function(gr,bpPastTSS=2000){
  # filter out genes that have no remainder
  gr <- gr[width(gr)>bpPastTSS]

  # filter out genes that have no strand (undefined behavior)
  gr <- gr[strand(gr)!="*"]

  #deal with genes on positive strand
  gr[strand(gr)=="+"] <- narrow(gr[strand(gr)=="+"],start=(bpPastTSS+1))

  #deal with genes on negative strand
  gr[strand(gr)=="-"] <- narrow(gr[strand(gr)=="-"],end=-1*(bpPastTSS+1))

  gr
} # end removePromRegion function

# function to add normalized square root of sum of squares of track to grange
addSqrTrackToGR <- function(inputGR,GRpath,GRfilename,fileLabel,GUsize=150){
  if(substr(GRfilename,nchar(GRfilename),nchar(GRfilename))%in% c("w","g","W","G")){
    GR <- import.bw(file.path(GRpath,GRfilename,fsep=""))
  }else{
    GR <- import.bedGraph(file.path(GRpath,GRfilename,fsep=""))
  }
  genome(GR) <- genome(inputGR)
  GR <- sortSeqlevels(GR)
  GR$score <- ((GR$score)^2)/GUsize

  score <- coverage(GR,weight="score")
  scoreCov <- coverage(GR,weight=(1./GUsize))
  inputGR <- binnedSumms(inputGR,score,paste("score_",fileLabel,sep=""))
  inputGR <- binnedSumms(inputGR,scoreCov,paste("numGUs_",fileLabel,sep=""))

  mcols(inputGR)[[paste("normedScore_",fileLabel,sep="")]] <- sqrt(
    mcols(inputGR)[[paste("score_",fileLabel,sep="")]] /
      mcols(inputGR)[[paste("numGUs_",fileLabel,sep="")]] )
  mcols(inputGR)[[paste("score_",fileLabel,sep="")]] <- NULL

  inputGR
} # end addSqrTrackToGR function

# helper function to find index of gene that has matching length
indexOfLengthMatch <- function(x,logLength){
  which(abs(x-logLength)==min(abs(x-logLength)))[1]
}

# Function to compute p-values using gammalss package.
getPvaluesRegions <- function(refVrefFiles,testVrefFiles,geneRegions,regionTitle,
                                  inFolder,outFolder,tName,rName,minGUs=10){

  #####################################################
  # set up geneRegions with data
  #####################################################

  # load null tracks
  for (ind in 1:length(refVrefFiles)){
    geneRegions <- addSqrTrackToGR(geneRegions,inFolder,refVrefFiles[ind],paste(regionTitle,"_null",ind,sep=""))
  }

  for (ind in 1:length(testVrefFiles)){
    geneRegions <- addSqrTrackToGR(geneRegions,inFolder,testVrefFiles[ind],paste(regionTitle,"_alt",ind,sep=""))
  }

  #####################################################
  # make x and y vectors for regression modeling
  #####################################################
  xVec <- c()
  yVec <- c()
  for (ind in 1:length(refVrefFiles)){
    tempY <- mcols(geneRegions)[[paste("normedScore_",regionTitle,"_null",ind,sep="")]]
    tempNumGUs <- mcols(geneRegions)[[paste("numGUs_",regionTitle,"_null",ind,sep="")]]

    # find valid indices
    hasDataInds <- (tempNumGUs>=minGUs)
    isFiniteInds <- is.finite(tempY)
    isntMissingInds <- !is.na(tempY)
    isValidInds <- (hasDataInds&isFiniteInds&isntMissingInds)

    # append values with sufficient data to xVec and yVec
    xVec <- c(xVec,log2(tempNumGUs[isValidInds]))
    yVec <- c(yVec,tempY[isValidInds])

  }# end loop over null files in construction of xVec/yVec

  # arrage GRange data into dataframe for modeling
  df <- data.frame(y=yVec,x=xVec)
    #cat("Feature: ", regionTitle, ".\n")
    #cat("Reference replicate name: ", rName, ".\n")
    #cat("Number of pairs: ", nrow(df), ".\n")

  #####################################################
  # do logit SST modeling
  #####################################################

  # Add/Substract epsilon in case there are 0 or 1's
  df$y[df$y<=0] <- .Machine$double.eps
  df$y[df$y>=1] <- 1 - .Machine$double.eps

  if((max(df$x) - min(df$x) < 1) || (length(df$x) < 1000)){
    # Fit logitSST model without regressing (need to check!)
    logitSST_fit <- gamlss(y ~ 1, data = df, sigma.formula =~ 1,
                           nu.formula =~ 1, tau.formula =~ 1,
                           family = "logitSST", trace = TRUE, n.cyc = 30)
  }else{
    # Fit logitSST model regressing on gene size
    logitSST_fit <- gamlss(y ~ pb(x), data = df, sigma.formula = ~pb(x),
                           nu.formula = ~pb(x), tau.formula = ~pb(x),
                           family = "logitSST", trace = TRUE, n.cyc = 30)
  }

  #####################################################
  # make plots
  #####################################################

  # save plot of the fit of the regression
  pdf(file=paste(outFolder,rName,"-",regionTitle,"-logitSST-fit-residuals.pdf",sep=""))
  plot(logitSST_fit)
  dev.off()

  # save plot of the centiles
  pdf(file=paste(outFolder,rName,"-",regionTitle,"-logitSST-centiles.pdf",sep=""))
  centiles(logitSST_fit, xvar = df$x, xlab = expression(log[2](s)), ylab = 't')
  dev.off()

  #####################################################
  # for each null comparison, find the p-values
  #####################################################
  Pvals.hist <- c()
  for (ind in 1:length(refVrefFiles)){

    # add new pvalue column
    mcols(geneRegions)[[paste("pVal_",regionTitle,"_null",ind,sep="")]] <- 1.0

    # loop over all genes to find pVals
    for (geneInd in 1:length(geneRegions)){

      normedScoreRegion <- mcols(geneRegions)[[paste("normedScore_",regionTitle,"_null",ind,sep="")]][geneInd]
      numGUsRegion      <- mcols(geneRegions)[[paste("numGUs_",regionTitle,"_null",ind,sep="")]][geneInd]

      # check if data is sufficient
      if ((numGUsRegion >= minGUs )&(!is.na(normedScoreRegion))&(is.finite(normedScoreRegion))){

        # find the fitted index for the length of current gene
        indexWid <- indexOfLengthMatch(df$x,log2(numGUsRegion))

        # compute the p-value using the fitted model
        tempPval <- 1 - plogitSST(normedScoreRegion,
                                  mu = fitted(logitSST_fit,"mu")[indexWid],
                                  sigma = fitted(logitSST_fit,"sigma")[indexWid],
                                  nu = fitted(logitSST_fit,"nu")[indexWid],
                                  tau = fitted(logitSST_fit,"tau")[indexWid])
        mcols(geneRegions)[[paste("pVal_",regionTitle,"_null",ind,sep="")]][geneInd] <- max(tempPval,.Machine$double.eps)
      }else{ # Not enough data, set to NA
        mcols(geneRegions)[[paste("pVal_",regionTitle,"_null",ind,sep="")]][geneInd] <- NA
      }# end check if enough data in geneBody

    }# end loop over geneRegions
    Pvals.hist <- c(Pvals.hist, mcols(geneRegions)[[paste("pVal_",regionTitle,"_null",ind,sep="")]])
  }# end loop over null files

  # Plot histogram of null p-values
  pdf(file=paste(outFolder,rName,"-",regionTitle,"-Pval-Null-Histogram.pdf",sep=""))
  hist(Pvals.hist, breaks = 100)
  dev.off()

  #####################################################
  # for each alternative comparison, find the p-values
  #####################################################
  Pvals.hist <- c()
  for (ind in 1:length(testVrefFiles)){

    # add new pvalue column
    mcols(geneRegions)[[paste("pVal_",regionTitle,"_alt",ind,sep="")]] <- 1.0

    # loop over all regions to find pVals
    for (geneInd in 1:length(geneRegions)){

      normedScoreRegion <- mcols(geneRegions)[[paste("normedScore_",regionTitle,"_alt",ind,sep="")]][geneInd]
      numGUsRegion     <- mcols(geneRegions)[[paste("numGUs_",regionTitle,"_alt",ind,sep="")]][geneInd]

      # if data is sufficient compute p-values
      if ((numGUsRegion >= minGUs )&(!is.na(normedScoreRegion))&(is.finite(normedScoreRegion))){

        # find the fitted index for the length of current gene
        indexWid <- indexOfLengthMatch(df$x,log2(numGUsRegion))

        # compute the p-value using the fitted model
        tempPval <-  1 - plogitSST(normedScoreRegion,
                                   mu = fitted(logitSST_fit,"mu")[indexWid],
                                   sigma =fitted(logitSST_fit,"sigma")[indexWid],
                                   nu = fitted(logitSST_fit,"nu")[indexWid],
                                   tau= fitted(logitSST_fit,"tau")[indexWid])

        mcols(geneRegions)[[paste("pVal_",regionTitle,"_alt",ind,sep="")]][geneInd] <- max(tempPval,.Machine$double.eps)

      }else{ # Not enough data, set to NA
        mcols(geneRegions)[[paste("pVal_",regionTitle,"_alt",ind,sep="")]][geneInd] <- NA
      }# end check if enough data in geneBody

    }# end loop over geneRegions
    Pvals.hist <- c(Pvals.hist, mcols(geneRegions)[[paste("pVal_",regionTitle,"_alt",ind,sep="")]])
  }# end loop over alternative files

  # Plot histogram of alt p-values
  pdf(file=paste(outFolder,rName,"-VS-",tName,"-",regionTitle, "Pval-Alt-Histogram.pdf", sep = ""))
  hist(Pvals.hist, breaks = 100)
  dev.off()

  #####################################################
  # return a list with the computed values
  #####################################################
  list(model = logitSST_fit, nullData = df, geneRegions = geneRegions)

} # end getPvaluesRegions function

rankBreakTies <- function(firstCol, secCol)
{
  # Return vector initizalized with string NA so that it doesn't become ND later on
  return_rank <- vector(length = length(firstCol))
  return_rank[1:length(return_rank)] <- "NA"
  
  # Ranking based on first column with min option
  rank1 <- BiocGenerics::rank(round(firstCol,digits = 16), ties.method = "min")
  
  # Get unique rankings
  rank1_unique <- BiocGenerics::unique(rank1)
   
  # For each unique rank, break ties with second column if data in first column
  N_data <- length(which(!is.na(firstCol)))
  for (r in rank1_unique) {
    # If ranking is greater than # of points with data then leave it as NA
    if(r>N_data) next
    # Subset data
    index <- which(rank1 == r)
    firstCol_sub <- firstCol[index]
    secCol_sub <- secCol[index]
    # If only more than one with that ranking, break tie
    if(length(firstCol_sub>1)){
      # Ranking offset given by second column - 1
      rank2 <- BiocGenerics::rank(round(secCol_sub,digits = 16), ties.method = "min") - 
        1 + rank1[index]
      return_rank[index] <- rank2
    } else {
      # No need to break tie
      return_rank[index] <- r
    }
  }
  # Return ranking
  return(return_rank)
} # end rankBreakTies function

rankGeneBodiesAndProms <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,
                                   tName,rName,chrsOfInterest = paste("chr",1:22,
                                   sep=""),bpPastTSS=2000,bpBeforeTSS=2000){

  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }

  ###############################################################
  # compute genes and  promoters
  ###############################################################
  # Get genes with entrez ID annotation
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb)
  genes <- genes[seqnames(genes) %in% chrsOfInterest]
  seqlevels(genes) <- seqlevels(genes)[1:length(chrsOfInterest)]

  # compute promoters
  proms <- promoters(genes,upstream=bpBeforeTSS,downstream=bpPastTSS)
  proms <- proms[seqnames(proms) %in% chrsOfInterest]
  seqlevels(proms) <- seqlevels(proms)[1:length(chrsOfInterest)]

  # run function to truncate genebodies to exclude promoters
  genes <- removePromRegion(genes,bpPastTSS=bpPastTSS)

  ###############################################################
  #run function to analyze promoters
  ###############################################################
  promsReturnObj <- getPvaluesRegions(refVrefFiles,testVrefFiles,proms,"prom",
                                      inFolder,outFolder,tName,rName)

  ###############################################################
  # run function to analyze truncated genebodies
  ###############################################################
  bodiesReturnObj <- getPvaluesRegions(refVrefFiles,testVrefFiles,genes,"bod",
                                       inFolder,outFolder,tName,rName)

  ###############################################################
  # combine promoter and bodies DF
  ###############################################################
  promsDF  <- as.data.frame(promsReturnObj$geneRegions)
  bodiesDF <- as.data.frame(bodiesReturnObj$geneRegions)

  #combine into a single DF with pVals
  mergedDF <- merge(promsDF[,c(grep("gene_id",names(promsDF)),
                    grep("^pVal", names(promsDF) ))],
                    bodiesDF[,c(grep("gene_id",names(bodiesDF)),
                    grep("^pVal", names(bodiesDF) ))],
                    by="gene_id",all=TRUE)
    # library(dplyr)
    # mergedDF <- full_join(promsDF[,c(grep("gene_id",names(promsDF)),
    #                       grep("^pVal", names(promsDF) ))],
    #                       bodiesDF[,c(grep("gene_id",names(bodiesDF)),
    #                       grep("^pVal", names(bodiesDF) ))],
    #                       by="gene_id")

  rm(list=c("promsDF","bodiesDF","promsReturnObj","bodiesReturnObj"))

  ###############################################################
  # compute combined prom/body statistics Tpb
  ###############################################################
  for (ind in 1:length(refVrefFiles)){
    mergedDF[[paste("combStat_null",ind,sep="")]] <-
      - 2 * log(mergedDF[,paste("pVal_prom_null",ind,sep="")]) -
      2 * log(mergedDF[,paste("pVal_bod_null",ind,sep="")])
  }
  for (ind in 1:length(testVrefFiles)){
    mergedDF[[paste("combStat_alt",ind,sep="")]] <-
      - 2 * log(mergedDF[,paste("pVal_prom_alt",ind,sep="")]) -
      2 * log(mergedDF[,paste("pVal_bod_alt",ind,sep="")])
  }

  ###############################################################
  # P-value, Q-value, and ranking of promoters
  ###############################################################
  for (ind in 1:nrow(mergedDF)){
    pVals <- mergedDF[ind,paste("pVal_prom_alt",1:length(testVrefFiles),sep="")]
    # do fisher's if there is data
    numNotNA <- sum(!is.na(pVals))
    if (numNotNA>1){
      mergedDF[ind,"pVal_prom_total"] <- pchisq(-2*sum(log(pVals),na.rm=TRUE),
        2*numNotNA, lower.tail=FALSE)
    } else if(numNotNA==1){
      mergedDF[ind,"pVal_prom_total"] <- pVals[!is.na(pVals)]
    } else{
      mergedDF[ind,"pVal_prom_total"] <- NA
    }
  }
  mergedDF$pVal_prom_total[which(mergedDF$pVal_prom_total < .Machine$double.eps)] <- .Machine$double.eps
  mergedDF$BH_qValue_prom <- p.adjust(mergedDF$pVal_prom_total,method="BH")

  ###############################################################
  # P-value, Q-value, and ranking of gene bodies
  ###############################################################
  for (ind in 1:nrow(mergedDF)){
    pVals <- mergedDF[ind,paste("pVal_bod_alt",1:length(testVrefFiles),sep="")]
    # do fisher's if there is data
    numNotNA <- sum(!is.na(pVals))
    if (numNotNA>1){
      mergedDF[ind,"pVal_bod_total"] <- pchisq( -2*sum(log(pVals),na.rm=TRUE),
        2*numNotNA, lower.tail=FALSE)
    } else if(numNotNA==1){
      mergedDF[ind,"pVal_bod_total"] <- pVals[!is.na(pVals)]
    } else{
      mergedDF[ind,"pVal_bod_total"] <- NA
    }
  }
  mergedDF$pVal_bod_total[which(mergedDF$pVal_bod_total < .Machine$double.eps)] <- .Machine$double.eps
  mergedDF$BH_qValue_bod <- p.adjust(mergedDF$pVal_bod_total,method="BH")

  #####################################################
  # make x and y vectors for null of Tpb
  #####################################################
  nullVec <- c()
  for (ind in 1:length(refVrefFiles)){
    tempY <- mergedDF[[paste("combStat_null",ind,sep="")]]
    notNAind <- !is.na(tempY)
    notINFind <- !is.infinite(tempY)

    # append values to nullVec
    nullVec <- c(nullVec, tempY[notNAind&notINFind])
  }# end loop over null files in construction of xVec/nullVec

  # compute null dist of prom/body stats
  pVals <- constructNullPvals(nullVec)

  #####################################################
  # make plots
  #####################################################
  # save plot of the empirical distribution versus theoretical chi-square
  pdf(file=paste(outFolder, tName, "-VS-", rName,
    "-emp-VS-chi2-bodyPromCombined.pdf", sep = ""))
  hist(nullVec, breaks = 250, prob = TRUE)
  curve(dchisq(x, df = 4), col = 'green', add = TRUE)
  dev.off()

  #####################################################
  # for each alt comparison, find the p-values
  #####################################################
  for (ind in 1:length(testVrefFiles)){

    # add new pvalue column
    mergedDF[[paste("pVal_comb_alt",ind,sep="")]] <- 1.0

    # deal with missing data in promoter (i.e., take value from body)
    missingProms <- is.na(mergedDF[[paste("pVal_prom_alt",ind,sep="")]])
    mergedDF[missingProms,paste("pVal_comb_alt",ind,sep="")] <-
      mergedDF[missingProms,paste("pVal_bod_alt",ind,sep="")]

    # deal with missing data in body (i.e., take value from promoter)
    missingBods <- is.na(mergedDF[[paste("pVal_bod_alt",ind,sep="")]])
    mergedDF[missingBods,paste("pVal_comb_alt",ind,sep="")] <-
      mergedDF[missingBods,paste("pVal_prom_alt",ind,sep="")]

    # compute pVals over remaining genes
    genesWithData <-  !(missingBods|missingProms)
    mergedDF[[paste("pVal_comb_alt",ind,sep="")]][genesWithData] <-
      pVals( mergedDF[[paste("combStat_alt",ind,sep="")]][genesWithData] )

  }

  ###############################################################
  # combine p-values using Fisher's meta-analysis
  ###############################################################
  for (ind in 1:nrow(mergedDF)){
    pVals <- mergedDF[ind,paste("pVal_comb_alt",1:length(testVrefFiles),sep="")]

    # do fisher's if there is data
    numNotNA <- sum(!is.na(pVals))
    if (numNotNA>1){
      mergedDF[ind,"pVal_comb_total"] <- pchisq( -2*sum(log(pVals),na.rm=TRUE),
        2*numNotNA, lower.tail=FALSE)
    } else if(numNotNA==1){
      mergedDF[ind,"pVal_comb_total"] <- pVals[!is.na(pVals)]
    } else{
      mergedDF[ind,"pVal_comb_total"] <- NA
    }
  }

  ###############################################################
  # return computed values sorted by final pVal
  ###############################################################
  rownames(mergedDF) <- mergedDF$gene_id
  ord <- order(mergedDF$pVal_comb_total)

  # set minimum p-value to be Machine precision
  mergedDF$pVal_comb_total[which(mergedDF$pVal_comb_total < .Machine$double.eps)] <- .Machine$double.eps

  # sort table
  mergedDF <- mergedDF[ord,]

  # find names of genes
  nameMapping <- AnnotationDbi::select(Homo.sapiens,key=as.character(mergedDF$gene_id),
                                       keytype="ENTREZID",columns=c("SYMBOL"))
  mergedDF$gene_name <- nameMapping$SYMBOL

  # Filter out repeated genes or NAs
  mergedDF <- mergedDF[!is.na(mergedDF$gene_name),]
  mergedDF <- mergedDF[!duplicated(mergedDF$gene_name),]
  mergedDF$BH_qValue <- p.adjust(mergedDF$pVal_comb_total,method="BH")

  ###############################################################
  # break any ties using rank prod
  ###############################################################
  mergedDF$rankProd    <- 1.0
  mergedDF$rankProdBod <- 1.0
  mergedDF$rankProdProm <- 1.0
  for (ind in 1:length(testVrefFiles)){
    ranks <- rank(mergedDF[[paste("pVal_comb_alt",ind,sep="")]])/
      length(mergedDF[[paste("pVal_comb_alt",ind,sep="")]])
    mergedDF$rankProd <- (mergedDF$rankProd)*ranks

    ranks <- rank(mergedDF[[paste("pVal_bod_alt",ind,sep="")]])/
      length(mergedDF[[paste("pVal_bod_alt",ind,sep="")]])
    mergedDF$rankProdBod <- (mergedDF$rankProdBod)*ranks

    ranks <- rank(mergedDF[[paste("pVal_prom_alt",ind,sep="")]])/
      length(mergedDF[[paste("pVal_prom_alt",ind,sep="")]])
    mergedDF$rankProdProm <- (mergedDF$rankProdProm)*ranks
  }

  # Ranking of promoters and gene bodies using respective P's first
  mergedDF$promRanks <- rankBreakTies(as.numeric(mergedDF$pVal_prom_total),
                                      as.numeric(mergedDF$rankProdProm))
  mergedDF$bodRanks <- rankBreakTies(as.numeric(mergedDF$pVal_bod_total),
                                     as.numeric(mergedDF$rankProdBod))

  # sort table
  mergedDF$Rank <- rankBreakTies(as.numeric(mergedDF$pVal_comb_total),
                                 as.numeric(mergedDF$rankProd))
  mergedDF <- mergedDF[order(mergedDF$Rank),]
  
  # Mount output table
  outTable <- cbind(mergedDF$gene_name,
                    mergedDF$Rank,
                    mergedDF$pVal_comb_total,
                    mergedDF$BH_qValue,
                    mergedDF[,c(grep("^pVal_comb_alt.*",names(mergedDF)))],
                    mergedDF$promRanks,
                    mergedDF$pVal_prom_total,
                    mergedDF$BH_qValue_prom,
                    mergedDF[,c(grep("^pVal_prom_alt.*",names(mergedDF)))],
                    mergedDF$bodRanks,
                    mergedDF$pVal_bod_total,
                    mergedDF$BH_qValue_bod,
                    mergedDF[,c(grep("^pVal_bod_alt.*",names(mergedDF)))]
                    )
  colnames(outTable) <- c("Gene",
                          "RANK comp.",
                          "p-value comp.",
                          "q-value comp.",
                          paste("p-value comp. TR",seq(1,length(testVrefFiles)),sep=""),
                          "RANK prom.",
                          "p-value prom.",
                          "q-value prom.",
                          paste("p-value prom. TR",seq(1,length(testVrefFiles)),sep=""),
                          "Rank_body",
                          "p-value body",
                          "q-value body",
                          paste("p-value body TR",seq(1,length(testVrefFiles)),sep="")
                          )
  rm(list=c("mergedDF"))
  
  # Remove NAs
  outTable[is.na(outTable)]<-"ND"
  
  ###############################################################
  # write results to output
  ###############################################################
  write.table(outTable,file=paste(outFolder,"gRank-",tName,"-VS-",rName,
              ".tsv",sep=""),row.names = FALSE,quote = FALSE,sep="\t")

  # return final dataframe
  outTable

} # end rankGeneBodiesAndProms function

# Define function that takes annotation file (BED,GFF,bedGraph,WIG,BigWig)
# and computes p-values. Builds null model, computes p-values, does Fisher's
# meta-analysis and rank product to rank the annotated features.
rankRegions <- function(refVrefFiles,testVrefFiles,regionsFile,regionsName,inFolder,
                        outFolder,tName="test",rName="ref",minGUs=10,
                        chrsOfInterest=paste("chr",1:22,sep=""))
{

  ###############################################################
  # check IO
  ###############################################################

  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }

  ###############################################################
  # retrieve seqinfo
  ###############################################################

  # Retrieve information about genome from UCSC
  genome.seqinfo <- Seqinfo(genome="hg19")
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  regions.GR <- rtracklayer::import(file.path(inFolder,regionsFile,fsep=""))
  regions.GR <- regions.GR[seqnames(regions.GR) %in% chrsOfInterest]
  genome(regions.GR) <- genome(genome.seqinfo)
  seqlevels(regions.GR) <- seqlevels(genome.seqinfo)
  seqlengths(regions.GR) <- seqlengths(genome.seqinfo)

  ###############################################################
  # run function to analyze regions
  ###############################################################

  # Compute p-values and generate GoF plots
  returnObj <- getPvaluesRegions(refVrefFiles,testVrefFiles,regions.GR,regionsName,
                                 inFolder,outFolder,tName,rName,minGUs)

      # Debugging
      # returnObj.backup <- returnObj

  ###############################################################
  # annotate with closest gene
  ###############################################################

  # Retrieve genes from UCSC
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  chr.filter <- list(tx_chrom = chrsOfInterest)
  annots.genes.gr <- genes(txdb,filter=chr.filter)
  symbol.table <- select(Homo.sapiens,key=annots.genes.gr$gene_id,keytype="ENTREZID",columns=c("SYMBOL"))
  annots.genes.gr$symbol <- symbol.table$SYMBOL
  annots.genes.gr <- unique(annots.genes.gr[!is.na(annots.genes.gr$symbol),])

  # Find closest gene
  returnObj.GR <- returnObj$geneRegions
  olaps.genes.gr <- IRanges::nearest(returnObj.GR,annots.genes.gr,select="arbitrary",ignore.strand=TRUE)
  returnObj.GR$symbol <- annots.genes.gr[olaps.genes.gr,]$symbol
  returnObj.GR$distance <- distance(returnObj.GR,annots.genes.gr[olaps.genes.gr,])
  returnObj.DF <- as.data.frame(returnObj.GR)

     # Debugging
     # returnObj.DF.backup <- returnObj.DF

  ###############################################################
  # combine p-values using Fisher's meta-analysis
  ###############################################################

  pval.col.alt.prefix <- paste("pVal_",regionsName,"_alt",sep="")
  pval.col.tot.prefix <- paste("pVal_",regionsName,"_total",sep="")
  for (ind in 1:nrow(returnObj.DF)){
    pVals <- returnObj.DF[ind,paste(pval.col.alt.prefix,1:length(testVrefFiles),sep="")]
    numNotNA <- sum(!is.na(pVals))
    if (numNotNA>1){
      returnObj.DF[ind,pval.col.tot.prefix] <- pchisq( -2*sum(log(pVals),na.rm=TRUE), 2*numNotNA, lower.tail=FALSE)
    } else if(numNotNA==1){
      returnObj.DF[ind,pval.col.tot.prefix] <- pVals[!is.na(pVals)]
    } else{
      returnObj.DF[ind,pval.col.tot.prefix] <- NA
    }
  }

  # Correct for MH
  returnObj.DF[which(returnObj.DF[,pval.col.tot.prefix] < .Machine$double.eps), pval.col.tot.prefix] <- .Machine$double.eps
  returnObj.DF$BH_qValue <- p.adjust(returnObj.DF[,pval.col.tot.prefix],method="BH")
  
  ###############################################################
  # break any ties using rank prod
  ###############################################################
  
  # Break any ties using rank prod
  returnObj.DF$rankProd    <- 1.0
  for (ind in 1:length(testVrefFiles)){
    ranks <- rank(returnObj.DF[[paste(pval.col.alt.prefix,ind,sep="")]])/
      length(returnObj.DF[[paste(pval.col.alt.prefix,ind,sep="")]])
    returnObj.DF$rankProd <- (returnObj.DF$rankProd)*ranks
  }
  
  # Identify NAs (problem with 2 arguments in order)
  returnObj.DF$Rank <- rankBreakTies(as.numeric(returnObj.DF[,pval.col.tot.prefix]),
                                     as.numeric(returnObj.DF$rankProd))
  returnObj.DF <- returnObj.DF[order(returnObj.DF$Rank),]
  
  # Remove NAs
  returnObj.DF[is.na(returnObj.DF)]<-"ND"
  
  # Construct output table. Need to substract 1 from start to make it 1-based coordinates.
  outTable <- cbind(
    returnObj.DF$seqnames,
    returnObj.DF$start-1,
    returnObj.DF$end,
    returnObj.DF$symbol,
    returnObj.DF$distance,
    returnObj.DF$Rank,
    returnObj.DF[,pval.col.tot.prefix],
    returnObj.DF$BH_qValue,
    returnObj.DF[,c(grep(paste("pVal_",regionsName,"_alt",sep=""),names(returnObj.DF)))]
  )
  colnames(outTable) <- c(
    "chromosome",
    "start site",
    "end site",
    "nearest gene",
    "DIST",
    "RANK",
    "p-value",
    "q-value",
    paste("p-value TR",seq(1,length(testVrefFiles)),sep="")
  )
  
  ###############################################################
  # write results to output
  ###############################################################

  # Write table
  write.table(outTable,file=paste(outFolder,"rRank-",regionsName,"-",
              rName,"-VS-",tName,".tsv",sep=""),row.names = FALSE,
              quote = FALSE,sep="\t")

  # return final dataframe
  outTable

} # end rankRegions function
