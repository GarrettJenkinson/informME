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
# Last Modified: 02/10/2017
# *************************************************************************
#
# REQUIRED PACKAGES:
# rtracklayer
# logitnorm
# mixtools

# Define smoothing function.
doSmoothing <- function(file,inFolder,outFolder,chrsOfInterest=paste("chr",1:22,sep=""),bandwidthVal=50000,outflag=FALSE) {
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
    gr.chr$score <- ksmooth(x=start(gr.chr),
                            y=gr.chr$score,
                            kernel="normal",
                            bandwidth=bandwidthVal,
                            x.points=start(gr.chr))$y
    gr <- c(gr,gr.chr)
  }
  gr <- gr[!is.na(gr$score)]
  
  if (length(gr)>0 && outflag){ # To output, set outflag=TRUE.
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("s",file,sep=""), 
                       visibility="full",autoScale=FALSE,viewLimits=c(0,1.0))
    export.bedGraph(gr,file.path(outFolder,paste("s",file,sep=""),fsep = ""),
                    trackLine=myTrackLine)
  }
  
  gr
} # End smoothing.

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

# Define thresholding and morphological closing function.
# This function operates on gr$score. 
doThreshMorph <- function(gr,file,outFolder,threshVal=50,
                          bandwidthVal=50000,GUsize=150.0,
                          requiredPercentBand=0.5){
  suppressMessages(library(rtracklayer))
  
  # Check folders for trailing slash, add if missing.
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Remove non-available (na) and infinite values.
  gr <- gr[!is.na(gr$score)]
  gr$score[is.infinite(gr$score)] <- 5*threshVal
  
  # Remove ranges below threshVal.
  grThresh <- gr[gr$score>threshVal]
  
  if (length(grThresh)>0){
    #
    # Do morphological closing.
    #
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
    
    #
    # Sum values inside the grThresh objects.
    #
    score <- coverage(gr,weight="score")
    cov   <- coverage(gr)
    grThresh <- binnedSumms(grThresh,score,"score")
    grThresh$score <- as.numeric(grThresh$score)/GUsize
    grThresh <- binnedSumms(grThresh,cov,"coverage")
    
    # remove DMRs that do not have enough data and get rid of coverage column
    grThresh <- grThresh[grThresh$coverage >= (bandwidthVal*requiredPercentBand)]
    grThresh$coverage <- NULL  
    
    #
    # Write output to file.
    #
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("DMR-",file,sep=""), 
                       visibility="full",autoScale=TRUE)
    export.bedGraph(grThresh[!is.na(grThresh$score)],file.path(outFolder,paste("DMR-",file,sep=""),fsep = ""),
                    trackLine=myTrackLine)
  } else{
    print("No DMRs - no output file is generated for:")
    print(file.path(outFolder,paste("DMR-",file,sep=""),fsep=""))
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

# Perform MH testing: BH or BY
multipleHypothesis <- function(nullGRs,altGRs,numNullComp,numAltComp,correction){
  # Collapse all pvalues
  x <- c()
  for(ind in 1:numNullComp){
    x <- rbind(x,data.frame(group='null',index=ind,pVals=nullGRs[[ind]]$pVals))
  }
  for(ind in 1:numAltComp){
    x <- rbind(x,data.frame(group='alt',index=ind,pVals=altGRs[[ind]]$pVals))
  }
  
  # Adjust p-values
  if(correction=='BH'){
    # Adjust p-values based on BH
    x$qVals <- p.adjust(x$pVals, method = 'BH')
  } else if(correction=='BY'){
    # Adjust p-values based on BH
    x$qVals <- p.adjust(x$pVals, method = 'BY')
  } else {
    # Correction type not valid
    print("Multiple hypothesis correction not valid.")
    exit(1)   
  }
  
  # Correct q-values smaller than machine precision.
  x$qVals[x$qVals<.Machine$double.eps] <- .Machine$double.eps
  
  # Add q-values with original GR objects
  for(ind in 1:numNullComp){
    nullGRs[[ind]]$qVals <- x[(x$group=='null') & (x$index==ind),]$qVals
  }
  for(ind in 1:numAltComp){
    altGRs[[ind]]$qVals <- x[(x$group=='alt') & (x$index==ind),]$qVals
  }

  # Return
  list(nullGRs,altGRs) 
}

# Function to compute p-values from mixture of normals in logit space.
# Used when replicate refernce data is not available.
logitMixturePvals <- function(values){
  suppressMessages(library(logitnorm))
  suppressMessages(library(mixtools))
  
  # Save 0,1 boundary values for the end.
  maxVals <- values>=1
  minVals <- values<=0
  values[maxVals] <- NA
  values[minVals] <- NA
  
  # Do mixture modeling in logit space.
  maxiters <- 1000
  df <- data.frame(sJSDlogit=log((values)/(1-values)))
  mixmdl <- normalmixEM(df$sJSDlogit[!is.na(df$sJSDlogit)],mu=c(-2.0,0.0),sigma=c(0.5,0.5),maxit = maxiters)
  
  mu <- mixmdl$mu
  sigma <- mixmdl$sigma
  lambda <- mixmdl$lambda
  
  if (length(mixmdl$all.loglik)>=maxiters){ # Mixture modeling not converged.
                                            # Be on the safe side and use a 
                                            # conservative model with large 
                                            # mean/sigma.
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
  }else{ # Mixture modeling converged. 
         # Null model comes from mixture component with smaller mean.
    if (mu[1]<mu[2]){
      muNull <- mu[1]
      sigmaNull <- sigma[1]
    }else{
      muNull <- mu[2]
      sigmaNull <- sigma[2]
    }
  }
  
  # Compute p-values.
  pVals <- 1 - plogitnorm(values,mu=muNull,sigma=sigmaNull)
  pVals[maxVals] <- min(pVals,na.rm=TRUE) # 1 value maps to smallest pval.
  pVals[minVals] <- 1                     # 0 value maps to 1 pval.
  
  # Return.
  pVals
  
}# End logitMixturePvals.

# Function to perform DMR detection when replicate reference data 
# is available. 
runReplicateDMR <- function(refVrefFiles,testVrefFiles,inFolder,outFolder,FDR=0.05, 
                            maxSQS = 250, chrsOfInterest=paste("chr",1:22,sep=""), 
                            bandwidthVal=50000,correction='BH',outflag=FALSE) {

  # Check folders for trailing slash, add if missing
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Count number of null (reference/reference) comparisons and 
  # number of alternative (test/reference) comparisons.
  numNullComp <- length(refVrefFiles)
  numAltComp  <- length(testVrefFiles)
  
  # Initialize lists of genomic ranges.
  nullGRs <- list()
  altGRs  <- list()
  
  # Do smoothing.
  for(ind in 1:numNullComp){
    write(paste("[",date(),"]: Smoothing JSD reference sample",ind,"out of",numNullComp), stderr())
    nullGRs[[ind]] <- doSmoothing(refVrefFiles[ind],inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,outflag=outflag)
  }
  for(ind in 1:numAltComp){
    write(paste("[",date(),"]: Smoothing JSD test sample",ind,"out of",numAltComp), stderr())
    altGRs[[ind]] <- doSmoothing(testVrefFiles[ind],inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,outflag=outflag)
  }
  
  # Find p-value function.
  nullVals <- vector()
  for (ind in 1:numNullComp){
    nullVals <- c(nullVals,nullGRs[[ind]]$score)
  }
  pValFn <- constructNullPvals(nullVals)
  
  # Find p-values.
  for(ind in 1:numNullComp){
    write(paste("[",date(),"]: Computing p-values reference sample",ind,"out of",numNullComp), stderr())
    nullGRs[[ind]]$pVals <- pValFn(nullGRs[[ind]]$score)
    # Correct p-values smaller than machine precision.
    nullGRs[[ind]]$pVals[nullGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  for(ind in 1:numAltComp){
    write(paste("[",date(),"]: Computing p-values test sample",ind,"out of",numAltComp), stderr())
    altGRs[[ind]]$pVals <- pValFn(altGRs[[ind]]$score)
    # Correct p-values smaller than machine precision.
    altGRs[[ind]]$pVals[altGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  
  # Estimate q-values if required: BH or BY
  if(!is.na(correction)){
    write(paste("[",date(),"]: Computing q-values based on",correction), stderr())
    out <- multipleHypothesis(nullGRs,altGRs,numNullComp,numAltComp,correction)
    nullGRs <- out[[1]]
    altGRs <- out[[2]]
  }

  # Write SQS files and DMRs.
  nullGRthresh <- list()
  altGRthresh  <- list()
  for (ind in 1:numNullComp){
    # Find SQS.
    nullGRs[[ind]]$sJSD <- nullGRs[[ind]]$score
    if(is.na(correction)){
      nullGRs[[ind]]$score <- -10*log10(nullGRs[[ind]]$pVals)
    } else {
      nullGRs[[ind]]$score <- -10*log10(nullGRs[[ind]]$qVals)
    }
    nullGRs[[ind]]$score[nullGRs[[ind]]$score>maxSQS] <- maxSQS

    # Generate SQS track    
    if(outflag){
      myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("SQS-s",refVrefFiles[ind],sep=""), 
                       visibility="full",autoScale=TRUE)
      export.bedGraph(nullGRs[[ind]],file.path(outFolder,paste("SQS-s",refVrefFiles[ind],sep=""),fsep=""))
    }

    # Morphological closing
    write(paste("[",date(),"]: Morphological closing reference sample",ind,"out of",numNullComp), stderr())
    if(is.na(correction)){
      nullGRthresh[[ind]] <- doThreshMorph(nullGRs[[ind]],refVrefFiles[ind],outFolder,bandwidthVal=bandwidthVal)
    } else {
      fdrToSqs <- -10*log10(FDR)
      nullGRthresh[[ind]] <- doThreshMorph(nullGRs[[ind]],refVrefFiles[ind],threshVal=fdrToSqs,outFolder,bandwidthVal=bandwidthVal)
    }
  }
  for (ind in 1:numAltComp){
    # Find SQS.
    altGRs[[ind]]$sJSD <- altGRs[[ind]]$score
    if(is.na(correction)){
      altGRs[[ind]]$score <- -10*log10(altGRs[[ind]]$pVals)
    } else {
      altGRs[[ind]]$score <- -10*log10(altGRs[[ind]]$qVals)
    }
    altGRs[[ind]]$score[altGRs[[ind]]$score>maxSQS] <- maxSQS

    # Generate SQS track    
    if(outflag){
      myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("SQS-s",testVrefFiles[ind],sep=""), 
                       visibility="full",autoScale=TRUE)
      export.bedGraph(altGRs[[ind]],file.path(outFolder,paste("SQS-s",testVrefFiles[ind],sep=""),fsep = ""))
    }
    
    # Morphological closing
    write(paste("[",date(),"]: Morphological closing test sample",ind,"out of",numAltComp), stderr())
    if(is.na(correction)){
      altGRthresh[[ind]] <- doThreshMorph(altGRs[[ind]],testVrefFiles[ind],outFolder,bandwidthVal=bandwidthVal)
    } else {
      fdrToSqs <- -10*log10(FDR)
      altGRthresh[[ind]] <- doThreshMorph(altGRs[[ind]],testVrefFiles[ind],threshVal=fdrToSqs,outFolder,bandwidthVal=bandwidthVal)
    }
  }
  
  # Return DMRs in alternative comparison.
  altGRthresh
} 

# Function to perform DMR detection when no replicate reference data 
# is available. 
runNoReplicateDMR <- function(file,inFolder,outFolder, maxSQS = 250,
                              chrsOfInterest=paste("chr",1:22,sep=""),
                              bandwidthVal=50000, correction='BH', 
			      FDR=0.05, outflag=FALSE) {
  
  
  # Check folders for trailing slash, add if missing.
  if(substr(inFolder,nchar(inFolder),nchar(inFolder)) != .Platform$file.sep ){
    inFolder <- paste(inFolder,.Platform$file.sep,sep="")
  }
  if(substr(outFolder,nchar(outFolder),nchar(outFolder)) != .Platform$file.sep ){
    outFolder <- paste(outFolder,.Platform$file.sep,sep="")
  }
  
  # Smooth track.
  write(paste("[",date(),"]: Smoothing JSD values"), stderr())
  GR <- doSmoothing(file,inFolder,outFolder,chrsOfInterest=chrsOfInterest,bandwidthVal=bandwidthVal,outflag=outflag)
  GR$sJSD <- GR$score
  
  # Find p-values.
  write(paste("[",date(),"]: Computing p-values"), stderr())
  GR$PvalsMix <- logitMixturePvals(GR$sJSD)
  
  # Adjust p-values
  if(correction=='BH'){
    write(paste("[",date(),"]: Computing q-values based on BH"), stderr())
    # Adjust p-values based on BH
    x$qVals <- p.adjust(x$pVals, method = 'BH')
    GR$score <- -10*log10(GR$qVals)
  } else if(correction=='BY'){
    write(paste("[",date(),"]: Computing q-values based on BY"), stderr())
    # Adjust p-values based on BH
    x$qVals <- p.adjust(x$pVals, method = 'BY')
    GR$score <- -10*log10(GR$qVals)
  } else if(is.na(correction)) {
    # Compute SQS score
    GR$score <- -10*log10(GR$PvalsMix)
  } else {
    # Correction type not valid
    print("Multiple hypothesis correction not valid.")
    exit(1)   
  }

  # Saturate SQS
  GR$score[GR$score>maxSQS] <- maxSQS
  
  # Save bedgraph file with SQS scores
  if(outflag){
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("SQS-s",file,sep=""), 
                       visibility="full",autoScale=TRUE)
    export.bedGraph(GR[!is.na(GR$score)],file.path(outFolder,paste("SQS-s",file,sep=""),fsep = ""))
  }
  
  # Do thresholding.
  write(paste("[",date(),"]: Morphological closing"), stderr())
  fdrToSqs <- -10*log10(FDR)
  GRthresh <- doThreshMorph(GR,file,outFolder,threshVal=fdrToSqs,bandwidthVal=bandwidthVal)
  
  # Return GR.
  GRthresh
}
