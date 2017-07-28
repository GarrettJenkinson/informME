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
# Last Modified: 07/28/2017
# *************************************************************************
#
# REQUIRED PACKAGES:
# rtracklayer
# logitnorm
# mixtools

#######################################################################################################################################
# GLOBAL VARIABLES
#######################################################################################################################################
PERM_CORR <- c('BH','BY')
SEED <- 100

#######################################################################################################################################
# FUNCTIONS
#######################################################################################################################################
# Define smoothing function.
doSmoothing <- function(file,inFolder,outFolder,chrsOfInterest=paste("chr",1:22,sep=""),bandwidthVal=50000,outflag=FALSE) {
  # Dependencies
  suppressMessages(library(rtracklayer))
  
  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Load data from file.
  GRTRACK <- import.bedGraph(file.path(inFolder,file,fsep=""),trackLine=FALSE)
  genome(GRTRACK) <- "hg19"
  GRTRACK <- sortSeqlevels(GRTRACK)
  
  # Go through chromosomes and smooth.
  gr <- GRTRACK[0]
  for (chr in chrsOfInterest){
    gr.chr <- GRTRACK[seqnames(GRTRACK) %in% chr]
    gr.chr <- sort(gr.chr)
    gr.chr$scoreOld <- gr.chr$score
    gr.chr$score <- ksmooth(x=start(gr.chr),y=gr.chr$score,kernel="normal",bandwidth=bandwidthVal,x.points=start(gr.chr))$y
    gr <- c(gr,gr.chr)
  }
  gr <- gr[!is.na(gr$score)]
  
  # To output, set outflag=TRUE.
  if (length(gr)>0 && outflag){ 
    myTrackLine<-new("GraphTrackLine",type="bedGraph",name=paste("s",file,sep=""),visibility="full",autoScale=FALSE,
		     viewLimits=c(0,1.0))
    export.bedGraph(gr,file.path(outFolder,paste("s",file,sep=""),fsep = ""),trackLine=myTrackLine)
  }
  
  # Return
  gr
} # End smoothing.

# Define function that computes the sum of SQS values (stored in numvar) within each DMR (stored in bins) - adapted from previously
# available documentation of the tileGenome function in GenomicRanges. 
binnedSumms <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

# Define thresholding and morphological closing function. This function operates on gr$score. 
doThreshMorph <- function(gr,file,outFolder,correction,pAdjThresh,bandwidthVal=50000,GUsize=150.0,requiredPercentBand=0.5){
  # Dependencies
  suppressMessages(library(rtracklayer))
  
  # Get equivalent SQS score
  sqsThreshVal <- -10*log10(pAdjThresh)

  # Check folders for trailing slash, add if missing.
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  basename <- paste("DMR-",correction,"-",pAdjThresh,"-",file,sep="")
  file.full.path <- file.path(outFolder,basename,fsep = "")
  
  # Remove non-available (na) and infinite values.
  gr <- gr[!is.na(gr$score)]
  gr$score[is.infinite(gr$score)] <- 5*sqsThreshVal
  
  # Remove ranges below sqsThreshVal.
  grThresh <- gr[gr$score>sqsThreshVal]
  
  if (length(grThresh[!is.na(grThresh$score)])>0){
    # Do morphological closing.
    strEl <- bandwidthVal
    
    # Increase size of intervals by 0.5*strEl on each side.
    grThreshUp <- flank(grThresh,round(0.5*strEl),start=TRUE)
    grThreshDown <- flank(grThresh,round(0.5*strEl),start=FALSE)
    
    # Compute union.
    grThresh <- reduce(c(grThresh,grThreshUp,grThreshDown))
    
    # Shrink size of intervals by 0.5*strEl on each side.
    grThreshUpSub <- flank(grThresh,round(-0.5*strEl),start=TRUE)
    grThreshDownSub <- flank(grThresh,round(-0.5*strEl),start=FALSE) 
    grThresh <- setdiff(grThresh,grThreshUpSub)
    grThresh <- setdiff(grThresh,grThreshDownSub)
    
    # Make sure everything is still within bounds.
    grThresh <- trim(grThresh)
    
    # Sum values inside the grThresh objects.
    score <- coverage(gr,weight="score")
    cov   <- coverage(gr)
    grThresh <- binnedSumms(grThresh,score,"score")
    grThresh$score <- as.numeric(grThresh$score)/GUsize
    grThresh <- binnedSumms(grThresh,cov,"coverage")
    
    # Remove DMRs that do not have enough data and get rid of coverage column
    grThresh <- grThresh[grThresh$coverage >= (bandwidthVal*requiredPercentBand)]
    grThresh$coverage <- NULL  
    
    if (length(grThresh[!is.na(grThresh$score)])>0){
      # Write output to file.
      myTrackLine <- new("GraphTrackLine",type="bedGraph",name=basename,visibility="full",autoScale=TRUE)
      export.bedGraph(grThresh[!is.na(grThresh$score)],file.full.path,trackLine=myTrackLine)
    } else {
      write(paste("[",date(),"]: No DMRs - no output file is generated for:"), stderr())
      write(paste("[",date(),"]: ",basename),stderr())
    }
  } else {
    write(paste("[",date(),"]: No DMRs - no output file is generated for:"), stderr())
    write(paste("[",date(),"]: ",basename),stderr())
  }
  
  # Return object.
  grThresh
  
}# End doThreshMorph function.

# Make empirical p-value function.
constructNullPvals <- function(nullValues){
  # Make null distribution.
  Fn <- ecdf(nullValues)
  
  # Return function to compute p-values.
  function(x) 1-Fn(x)
}

# Function to check correction method is valid
checkValidCorrection <- function(correction) {
  # Check if argument is in permitted set of corrections
  if(!(correction %in% PERM_CORR)){
    write(paste("[",date(),"]: Multiple hypothesis correction method",correction,"not currently supported."), stderr())
    write(paste("[",date(),"]: Must modify code to remove this error to proceed."), stderr())
    stop("Stopping due to invalid correction method")
  }
}

# Perform MH testing: Benjamini & Yekutieli (default) or Benjamini & Hochberg
multipleHypothesis <- function(nullGRs,altGRs,numNullComp,numAltComp,correction){
  # Collapse all pvalues
  x.null <- c()
  x.alt <- c()
  for(ind in 1:numNullComp){
    x.null <- rbind(x.null,data.frame(index=ind,pVals=nullGRs[[ind]]$pVals))
  }
  for(ind in 1:numAltComp){
    x.alt <- rbind(x.alt,data.frame(index=ind,pVals=altGRs[[ind]]$pVals))
  }
  
  # Adjust p-values
  x.null$pAdjVals <- p.adjust(x.null$pVals, method = correction)
  x.alt$pAdjVals <- p.adjust(x.alt$pVals, method = correction)

  # Correct q-values smaller than machine precision.
  x.null$pAdjVals[x.null$pAdjVals<.Machine$double.eps] <- .Machine$double.eps
  x.alt$pAdjVals[x.alt$pAdjVals<.Machine$double.eps] <- .Machine$double.eps
  
  # Add q-values with original GR objects
  for(ind in 1:numNullComp){
    nullGRs[[ind]]$pAdjVals <- x.null[x.null$index==ind,]$pAdjVals
  }
  for(ind in 1:numAltComp){
    altGRs[[ind]]$pAdjVals <- x.alt[x.alt$index==ind,]$pAdjVals
  }

  # Return
  list(nullGRs,altGRs) 
}

# Function to compute p-values from mixture of normals in logit space. Used when replicate refernce data is not available.
logitMixturePvals <- function(values){
  # Dependencies
  suppressMessages(library(logitnorm))
  suppressMessages(library(mixtools))
  
  # Save 0,1 boundary values for the end.
  maxVals <- values>=1
  minVals <- values<=0
  values[maxVals] <- NA
  values[minVals] <- NA
  
  # Do mixture modeling in logit space.
  write(paste("[",date(),"]: Running EM algorithm"), stdout())
  maxiters <- 1000
  df <- data.frame(sJSDlogit=log((values)/(1-values)))

  # EM Algorithm with fixed seed for reproducibility
  set.seed(SEED)
  mixmdl <- normalmixEM(df$sJSDlogit[!is.na(df$sJSDlogit)],mu=c(-2.0,0.0),sigma=c(0.5,0.5),maxit = maxiters)
 
  # Get parameters from mixture 
  mu <- mixmdl$mu
  sigma <- mixmdl$sigma
  lambda <- mixmdl$lambda
  
  # Check for convergence of the EM algorithm
  if (length(mixmdl$all.loglik)>=maxiters){ 
    # Mixture modeling not converged. Be on the safe side and use a conservative model with large mean/sigma.
    if (mu[1]>mu[2]){
      muNull <- mu[1]
    }else{
      muNull <- mu[2]
    }
    if (sigma[1]>sigma[2]){
      sigmaNull <- sigma[1]
    }else{
      sigmaNull <- sigma[2]
    }
  }else{ 
    # Mixture modeling converged. Null model comes from mixture component with smaller mean.
    if (mu[1]<mu[2]){
      muNull <- mu[1]
      sigmaNull <- sigma[1]
    }else{
      muNull <- mu[2]
      sigmaNull <- sigma[2]
    }
  }
  
  # Compute p-values.
  write(paste("[",date(),"]: Computing p-values"), stdout())
  pVals <- 1 - plogitnorm(values,mu=muNull,sigma=sigmaNull)
  pVals[maxVals] <- min(pVals,na.rm=TRUE) # 1 value maps to smallest pval.
  pVals[minVals] <- 1                     # 0 value maps to 1 pval.
  
  # Return.
  pVals
  
}# End logitMixturePvals.

# Function to perform DMR detection when replicate reference data is available. 
runReplicateDMR <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,maxSQS=250,chrsOfInterest=paste("chr",1:22,sep=""),
			    bandwidthVal=50000,correction='BY',pAdjThresh=0.01,outflag=FALSE) {

  # Check correction method
  checkValidCorrection(correction)

  # Check folders for trailing slash, add if missing
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Count number of null (reference/reference) comparisons and number of alternative (test/reference) comparisons.
  numNullComp <- length(refVrefFiles)
  numAltComp  <- length(testVrefFiles)
  
  # Initialize lists of genomic ranges.
  nullGRs <- list()
  altGRs  <- list()
  
  # Do smoothing.
  for(ind in 1:numNullComp){
    write(paste("[",date(),"]: Smoothing JSD reference sample",ind,"out of",numNullComp), stdout())
    nullGRs[[ind]] <- doSmoothing(refVrefFiles[ind],inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,
				  outflag=outflag)
  }
  for(ind in 1:numAltComp){
    write(paste("[",date(),"]: Smoothing JSD test sample",ind,"out of",numAltComp), stdout())
    altGRs[[ind]] <- doSmoothing(testVrefFiles[ind],inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,
			         outflag=outflag)
  }
  
  # Find p-value function.
  nullVals <- vector()
  for (ind in 1:numNullComp){
    nullVals <- c(nullVals,nullGRs[[ind]]$score)
  }
  pValFn <- constructNullPvals(nullVals)
  
  # Find p-values.
  for(ind in 1:numNullComp){
    write(paste("[",date(),"]: Computing p-values reference sample",ind,"out of",numNullComp), stdout())
    nullGRs[[ind]]$pVals <- pValFn(nullGRs[[ind]]$score)
    nullGRs[[ind]]$pVals[nullGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  for(ind in 1:numAltComp){
    write(paste("[",date(),"]: Computing p-values test sample",ind,"out of",numAltComp), stdout())
    altGRs[[ind]]$pVals <- pValFn(altGRs[[ind]]$score)
    altGRs[[ind]]$pVals[altGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  
  # Estimate q-values using Benjamini & Yekutieli or Benjamini & Hochberg. MH is done independently: one procedure for NULL, 
  # and a second one for ALT. 
  write(paste("[",date(),"]: Computing q-values based on",correction), stdout())
  out <- multipleHypothesis(nullGRs,altGRs,numNullComp,numAltComp,correction)
  nullGRs <- out[[1]]
  altGRs <- out[[2]]

  # Write SQS files and DMRs.
  nullGRthresh <- list()
  altGRthresh  <- list()
  
  # Null comparisons
  for (ind in 1:numNullComp){
    # Find SQS.
    nullGRs[[ind]]$sJSD <- nullGRs[[ind]]$score
    nullGRs[[ind]]$score <- -10*log10(nullGRs[[ind]]$pAdjVals)
    nullGRs[[ind]]$score[nullGRs[[ind]]$score>maxSQS] <- maxSQS

    # Generate SQS track    
    if(outflag){
      myTrackLine<-new("GraphTrackLine",type="bedGraph",name=paste("SQS-s",refVrefFiles[ind],sep=""),visibility="full",autoScale=TRUE)
      export.bedGraph(nullGRs[[ind]],file.path(outFolder,paste("SQS-s",refVrefFiles[ind],sep=""),fsep=""))
    }

    # Morphological closing
    write(paste("[",date(),"]: Morphological closing reference sample",ind,"out of",numNullComp), stdout())
    nullGRthresh[[ind]] <- doThreshMorph(nullGRs[[ind]],refVrefFiles[ind],outFolder,correction,pAdjThresh,bandwidthVal=bandwidthVal)
  }
  
  # Alternative comparisons
  for (ind in 1:numAltComp){
    # Find SQS.
    altGRs[[ind]]$sJSD <- altGRs[[ind]]$score
    altGRs[[ind]]$score <- -10*log10(altGRs[[ind]]$pAdjVals)
    altGRs[[ind]]$score[altGRs[[ind]]$score>maxSQS] <- maxSQS

    # Generate SQS track    
    if(outflag){
      myTrackLine<-new("GraphTrackLine",type="bedGraph",name=paste("SQS-s",testVrefFiles[ind],sep=""),visibility="full",autoScale=TRUE)
      export.bedGraph(altGRs[[ind]],file.path(outFolder,paste("SQS-s",testVrefFiles[ind],sep=""),fsep = ""))
    }
    
    # Morphological closing
    write(paste("[",date(),"]: Morphological closing test sample",ind,"out of",numAltComp), stdout())
    altGRthresh[[ind]] <- doThreshMorph(altGRs[[ind]],testVrefFiles[ind],outFolder,correction,pAdjThresh,bandwidthVal=bandwidthVal)
  }
  
  # Return DMRs in alternative comparison.
  altGRthresh
} 

# Function to perform DMR detection when no replicate reference data is available. 
runNoReplicateDMR <- function(file,inFolder,outFolder,maxSQS=250,chrsOfInterest=paste("chr",1:22,sep=""),bandwidthVal=50000,
			      correction='BY',pAdjThresh=0.01,outflag=FALSE) {
 
  # Check correction method
  checkValidCorrection(correction)
 
  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Smooth track.
  write(paste("[",date(),"]: Smoothing JSD values"), stdout())
  GR <- doSmoothing(file,inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,outflag=outflag)
  GR$sJSD <- GR$score
  
  # Find p-values.
  GR$PvalsMix <- logitMixturePvals(GR$sJSD)
  
  # Adjust p-values
  write(paste("[",date(),"]: Computing q-values based on ",correction," method."), stdout())
  GR$pAdjVals <- p.adjust(GR$PvalsMix, method = correction)
  GR$score <- -10*log10(GR$pAdjVals)

  # Saturate SQS
  GR$score[GR$score>maxSQS] <- maxSQS
  
  # Save bedgraph file with SQS scores
  if(outflag){
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("SQS-s",file,sep=""),visibility="full",autoScale=TRUE)
    export.bedGraph(GR[!is.na(GR$score)],file.path(outFolder,paste("SQS-s",file,sep=""),fsep = ""))
  }
  
  # Do thresholding.
  write(paste("[",date(),"]: Morphological closing"), stdout())
  GRthresh <- doThreshMorph(GR,file,outFolder,correction,pAdjThresh,bandwidthVal=bandwidthVal)
  
  # Return GR.
  GRthresh
}
