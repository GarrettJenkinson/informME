#
#   informME: An information-theoretic pipeline for WGBS data
#   Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software Foundation,
#   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#   or see <http://www.gnu.org/licenses/>.
#




#REQUIRED PACKAGES:
# rtracklayer
# logitnorm
# mixtools



# define function modified from ?tileGenome
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



# Define smoothing function:
doSmoothing <- function(file,inFolder,outFolder,chrsOfInterest=paste("chr",1:22,sep=""),bandwidthVal=50000) {
  library(rtracklayer)
  
  # load data from file
  GRTRACK <- import.bedGraph(paste(inFolder,file,sep=""),trackLine=FALSE)
  genome(GRTRACK) <- "hg19"
  GRTRACK <- sortSeqlevels(GRTRACK)
  
  # loop over chromosomes, smoothing each
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
  
  if (length(gr)>0){
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("smoothed-",file,sep=""), 
                       visibility="full",autoScale=FALSE,viewLimits=c(0,1.0))
    export.bedGraph(gr,paste(outFolder,"smoothed-",file,sep=""),
                    trackLine=myTrackLine)
  }
  
  gr
} # end doSmoothing function




# Define thresholding and morphological closing function.
# this function operates on gr$score 
doThreshMorph <- function(gr,file,outFolder,threshVal=50,bandwidthVal=50000,GUsize=150.0){
  library(rtracklayer)
  
  # remove na and infinite valeus
  gr <- gr[!is.na(gr$score)]
  gr$score[is.infinite(gr$score)] <- 5*threshVal
  
  # remove ranges below threshVal
  grThresh <- gr[gr$score>threshVal]
  
  if (length(grThresh)>0){
    #
    # Do morphological closing:
    #
    strEl <- (bandwidthVal-GUsize)
    
    # increase size of intervals by .5*strEl on each side
    grThreshUp <- flank(grThresh,round(.5*strEl),start=TRUE)
    grThreshDown <- flank(grThresh,round(.5*strEl),start=FALSE)
    
    # compute union
    grThresh <- reduce(c(grThresh,grThreshUp,grThreshDown))
    
    # shrink size of intervals by .5*strEl on each side
    grThreshUpSub <- flank(grThresh,round(-.5*strEl),start=TRUE)
    grThreshDownSub <- flank(grThresh,round(-.5*strEl),start=FALSE)
    
    grThresh <- setdiff(grThresh,grThreshUpSub)
    grThresh <- setdiff(grThresh,grThreshDownSub)
    
    # make sure everything is still in-bounds:
    grThresh <- trim(grThresh)
    
    #
    # Sum values inside the grThresh objects
    #
    score <- coverage(gr,weight="score")
    grThresh <- binnedSumms(grThresh,score,"score")
    grThresh$score <- as.numeric(grThresh$score)/GUsize
    
    #
    # write output to file
    #
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("thresh",threshVal,"-",file,sep=""), 
                       visibility="full",autoScale=TRUE)
    export.bedGraph(grThresh[!is.na(grThresh$score)],paste(outFolder,"thresh",threshVal,"-",file,sep=""),
                    trackLine=myTrackLine)
  } else{
    print("No DMRs, and so no file written for:")
    print(paste(outFolder,"thresh",threshVal,"-",file,sep=""))
  }
  
  #return object
  grThresh
  
}# end doThreshMorph function


#make the empirical p-value function
constructNullPvals <- function(nullValues){
  # make null dist
  Fn <- ecdf(nullValues)
  
  # return function to compute p-values
  function(x) 1-Fn(x)
}



# Function to compute p-values from mixture of normals in logit space.
# Used in the case when replicates are not available to construct a null.
logitMixturePvals <- function(values){
  library(logitnorm)
  library(mixtools)
  
  # save 0,1 boundary values for the end
  maxVals <- values>=1
  minVals <- values<=0
  values[maxVals] <- NA
  values[minVals] <- NA
  
  # do mixture modeling in logit space
  df <- data.frame(sJSDlogit=log((values)/(1-values)))
  mixmdl <- normalmixEM(df$sJSDlogit[!is.na(df$sJSDlogit)],mu=c(-2.0,0.0),sigma=c(0.5,0.5))
  
  mu <- mixmdl$mu
  sigma <- mixmdl$sigma
  lambda <- mixmdl$lambda
  
  if (mu[1]<mu[2]){
    muNull <- mu[1]
    sigmaNull <- sigma[1]
  }else{
    muNull <- mu[2]
    sigmaNull <- sigma[2]
  }
  
  # compute pvalues
  pVals <-  1 - plogitnorm(values,mu=muNull,sigma=sigmaNull)
  pVals[maxVals] <- min(pVals,na.rm=TRUE) # 1 value maps to smallest pval
  pVals[minVals] <- 1 # 0 value maps to 1 pval
  
  #return
  pVals
  
}# end logitMixturePvals


# function to run empirical null DMR algorithm
runEmpiricalNullDMR <- function(nullFiles,altFiles,inFolder,outFolder, maxSQS = 250) {
  
  # count the number of null (refeence/reference) comparisons and 
  # the number of alternative (test/reference) comparisons
  numNullComp <- length(nullFiles)
  numAltComp  <- length(altFiles)
  
  # initialize lists of genomic ranges
  nullGRs <- list()
  altGRs  <- list()
  
  # Do all smoothing
  for(ind in 1:numNullComp){
    nullGRs[[ind]] <- doSmoothing(nullFiles[ind],inFolder,outFolder)
  }
  for(ind in 1:numAltComp){
    altGRs[[ind]] <- doSmoothing(altFiles[ind],inFolder,outFolder)
  }
  
  # Find p-value function
  nullVals <- vector()
  for (ind in 1:numNullComp){
    nullVals <- c(nullVals,nullGRs[[ind]]$score)
  }
  pValFn <- constructNullPvals(nullVals)
  
  # Find p-values
  for(ind in 1:numNullComp){
    nullGRs[[ind]]$pVals <- pValFn(nullGRs[[ind]]$score)
    #correct pValues smaller than machine precision
    nullGRs[[ind]]$pVals[nullGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  for(ind in 1:numAltComp){
    altGRs[[ind]]$pVals <- pValFn(altGRs[[ind]]$score)
    #correct pValues smaller than machine precision
    altGRs[[ind]]$pVals[altGRs[[ind]]$pVals<.Machine$double.eps] <- .Machine$double.eps
  }
  
  # write out SQS files and DMRs
  nullGRthresh <- list()
  altGRthresh  <- list()
  for (ind in 1:numNullComp){
    #find SQS
    nullGRs[[ind]]$sJSD <- nullGRs[[ind]]$score
    nullGRs[[ind]]$score <-  -10*log10(nullGRs[[ind]]$pVals)
    nullGRs[[ind]]$score[nullGRs[[ind]]$score>maxSQS] <- maxSQS
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("smoothSQS-",nullFiles[ind],sep=""), 
                       visibility="full",autoScale=TRUE)
    export.bedGraph(nullGRs[[ind]],paste(outFolder,"smoothSQS-",nullFiles[ind],".bed",sep=""))
    
    nullGRthresh[[ind]] <- doThreshMorph(nullGRs[[ind]],paste("smoothSQS-",nullFiles[ind],sep=""),outFolder)
  }
  for (ind in 1:numAltComp){
    #find SQS
    altGRs[[ind]]$sJSD <- altGRs[[ind]]$score
    altGRs[[ind]]$score <-  -10*log10(altGRs[[ind]]$pVals)
    altGRs[[ind]]$score[altGRs[[ind]]$score>maxSQS] <- maxSQS
    myTrackLine <- new("GraphTrackLine",type="bedGraph", name=paste("smoothSQS-",altFiles[ind],sep=""), 
                       visibility="full",autoScale=TRUE)
    export.bedGraph(altGRs[[ind]],paste(outFolder,"smoothSQS-",altFiles[ind],".bed",sep=""))
    
    altGRthresh[[ind]] <- doThreshMorph(altGRs[[ind]],paste("smoothSQS-",altFiles[ind],sep=""),outFolder)
  }
  
  #return DMRs in alt comparison
  altGRthresh
} 


# function to run mixture logit null DMR algorithm
runMixLogitNullDMR <- function(file,inFolder,outFolder, maxSQS = 250) {
  
  # smooth the track
  GR <- doSmoothing(file,inFolder,outFolder)
  GR$sJSD <- GR$score
  
  # Find p-values
  GR$PvalsMix <- logitMixturePvals(GR$sJSD)
  GR$score <-  -10*log10(GR$PvalsMix)
  GR$score[GR$score>maxSQS] <- maxSQS
  myTrackLine  <-  new("GraphTrackLine",type="bedGraph", name=paste("smoothSQSmix-",file,sep=""), 
                       visibility="full",autoScale=TRUE)
  export.bedGraph(GR[!is.na(GR$score)],paste(outFolder,"smoothSQSmix-",file,sep=""))
  
  # do thresholding
  GRthresh <- doThreshMorph(GR,paste("SQSmix-",file7,sep=""),outFolder)
  
  #return GR
  GRthresh
}
