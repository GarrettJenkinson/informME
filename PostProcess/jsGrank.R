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
# *************************************************************************
# Last Modified: 04/17/2017
# *************************************************************************
#
# REQUIRED PACKAGES:
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("Homo.sapiens")

# Writing large excel files causes errors unless the Java parameters are 
# changed.
options(java.parameters = c("-Xmx6g","-XX:-UseConcMarkSweepGC") )
library("XLConnect")

# Define function to rank genes when no replicate reference data is available.
# Pass empty vector c() to refVrefFiles. 
rankGenes <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,
                      tName="test",rName="ref",
                      GUsize=150,upStreamSize=2000,downStreamSize=2000,
                      chrsOfInterest=paste("chr",1:22,sep=""),minGUs=1){
  
  if(length(refVrefFiles)==0){ # Use JSD to rank the genes.
    result <- rankGenesJSD(testVrefFiles,inFolder,outFolder,tName,rName,
                           GUsize=GUsize,upStreamSize=upStreamSize,downStreamSize=downStreamSize,
                           chrsOfInterest=chrsOfInterest,minGUs=minGUs)
    return(result)
  }
  # Else use empirical null distribution to rank the genes. 
  result <- rankGenesEmpNull(refVrefFiles,testVrefFiles,inFolder,outFolder,tName,rName,
                             GUsize=GUsize,upStreamSize=upStreamSize,downStreamSize=downStreamSize,
                             chrsOfInterest=chrsOfInterest,minGUs=minGUs)
}

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
}

# Define function to rank genes when replicate reference data is available. 
rankGenesEmpNull <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,tName="test",rName="ref",
                           GUsize=150,upStreamSize=2000,downStreamSize=2000,
                           chrsOfInterest=paste("chr",1:22,sep=""),minGUs=1,
                           doRankProd=FALSE){

  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Define promoter regions.
  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  genesAll<-transcripts(txdb)
  genes<-keepSeqlevels(genesAll,chrsOfInterest,pruning.mode="coarse")
  sortSeqlevels(genes)
  proms<-promoters(genes,upstream=upStreamSize,downstream=downStreamSize)
  
  # Count number of reference/reference comparisons and 
  # number of test/reference comparisons.
  numNullComp <- length(refVrefFiles)
  numAltComp  <- length(testVrefFiles)
  
  # Initialize lists of genomic ranges.
  nullGRs <- list()
  altGRs  <- list()
 
  # Load files into genomic range lists and arrange the seqlevels.  
  for(ind in 1:numNullComp){
    if(substr(refVrefFiles[ind],nchar(refVrefFiles[ind]),nchar(refVrefFiles[ind]))%in% c("w","g","W","G")){
      nullGRs[[ind]] <- import.bw(file.path(inFolder,refVrefFiles[ind],fsep=""))
    }else{
      nullGRs[[ind]] <- import.bedGraph(file.path(inFolder,refVrefFiles[ind],fsep=""))
    }
    genome(nullGRs[[ind]]) <- "hg19"
    sortSeqlevels(nullGRs[[ind]])
    nullGRs[[ind]] <- keepSeqlevels(nullGRs[[ind]],chrsOfInterest,pruning.mode="coarse")
  }
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
  
  # Build empirical null distribution.
  nullVals <- c()
  for (ind in 1:numNullComp){
    nullVals <- c(nullVals,nullGRs[[ind]]$score)
  }
  Fn <- ecdf(nullVals)
  
  # Find p-values and summand for Fisher's method.
  for(ind in 1:numAltComp){
    altGRs[[ind]]$pVals <- (1-Fn(altGRs[[ind]]$score)) # compute pValue
    
    # Correct for pValues smaller than a "reasonable" precision.
    altGRs[[ind]]$pVals[altGRs[[ind]]$pVals < .Machine$double.eps] <- .Machine$double.eps
    
    altGRs[[ind]]$score <- -2*log(altGRs[[ind]]$pVals) # Summand for Fisher's method.
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
    
    # Compute Fisher's meta analysis pvalue.
    mcols(proms)[[paste("metaPval",ind,sep="")]] <- pchisq(mcols(proms)[[paste("score",ind,sep="")]],
                                                           2*mcols(proms)[[paste("cov",ind,sep="")]],
                                                           lower.tail = FALSE)
    mcols(proms)[[paste("metaPval",ind,sep="")]][mcols(proms)[[paste("cov",ind,sep="")]]< minGUs] <- NA
  }
  
  # Generate ranked gene lists.
  rankingTabs <- list()
  for (ind in 1:numAltComp){
    sortedIndices <- order(mcols(proms)[[paste("metaPval",ind,sep="")]],decreasing = FALSE)
    sortedOutput <- select(Homo.sapiens, key=as.character(proms[sortedIndices]$tx_id), keytype="TXID", columns=c("SYMBOL"))
    sortedOutput <- sortedOutput[(!duplicated(sortedOutput$SYMBOL))&(!is.na(sortedOutput$SYMBOL)),,drop=FALSE]
    sortedOutput[[paste("rank",ind,sep="")]] <- 1:length(sortedOutput$SYMBOL)
    rankingTabs[[ind]] <- sortedOutput
    rm(list=c("sortedIndices","sortedOutput"))
  }
  rm(list="proms") 

  # Generate merged table. 
  df <- rankingTabs[[1]][,c(2,3)]
  if (numAltComp>1){
    for (ind in 2:numAltComp){
      df <- merge(df,rankingTabs[[ind]][,c(2,3)],all=TRUE)
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
  # Note that this requires the package rankprodbounds.R and an 
  # independence assumption that is not always valid.
  if (numAltComp>1 & doRankProd){
    # Get HesEisBre-14 algorithm for merging ranked lists.
    source(paste(inFolder,"rankprodbounds.R",sep=""))
    
    # Compute pValue.
    df$pVals<-rankprodbounds(df$rankprod,length(df$rankprod),numAltComp,Delta = 'geometric')
  }

  # Reformat table to handle one or more comparisons.  
  if (numAltComp>1){
    df$rankProd <- 1:length(df$rankProd)
  } else{
    df$rank <- df$rank1
    df <- subset(df,select = c(-rankProd,-rank1))  
  }
  
  # Write results to output.
  writeWorksheetToFile(paste(outFolder,"gRankRRD-JSD-",tName,"-VS-",rName,".xlsx",sep=""), 
                       data = df, sheet = "gRankRRD")
  gc() # Run garbage collection since dealing with xlsx is memory intensive. 
  
  # Return value.
  df
}# End rankGenesEmpNull.

# Define function to rank genes when no replicate reference data is available.
rankGenesJSD <- function(testVrefFiles,inFolder,outFolder,tName="test",rName="ref",
                             GUsize=150,upStreamSize=2000,downStreamSize=2000,
                             chrsOfInterest=paste("chr",1:22,sep=""),minGUs=1,
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
  genesAll<-transcripts(txdb)
  genes<-keepSeqlevels(genesAll,chrsOfInterest,pruning.mode="coarse")
  sortSeqlevels(genes)
  proms<-promoters(genes,upstream=upStreamSize,downstream=downStreamSize)
  
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
    sortedOutput <- select(Homo.sapiens, key=as.character(proms[sortedIndices]$tx_id), keytype="TXID", columns=c("SYMBOL"))
    sortedOutput <- sortedOutput[(!duplicated(sortedOutput$SYMBOL))&(!is.na(sortedOutput$SYMBOL)),,drop=FALSE]
    sortedOutput[[paste("rank",ind,sep="")]] <- 1:length(sortedOutput$SYMBOL)
    rankingTabs[[ind]] <- sortedOutput
    rm(list=c("sortedIndices","sortedOutput"))
  }
  rm(list="proms") 
  
  # Generate merged table. 
  df <- rankingTabs[[1]][,c(2,3)]
  if (numAltComp>1){
    for (ind in 2:numAltComp){
      df <- merge(df,rankingTabs[[ind]][,c(2,3)],all=TRUE)
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
    df$rankProd <- 1:length(df$rankProd)
  } else{
    df$rank <- df$rank1
    df <- subset(df,select = c(-rankProd,-rank1))  
  }
   
  # Write results to output.
  writeWorksheetToFile(paste(outFolder,"gRank-JSD-",tName,"-VS-",rName,".xlsx",sep=""), 
                       data = df, sheet = "gRank")
  gc() # Run garbage collection since dealing with xlsx is memory intensive. 
  
  # Return value.
  df
}# End rankGenesJSD.
