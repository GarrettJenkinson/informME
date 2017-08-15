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
  grThresh[!is.na(grThresh$score)]
  
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

# Function to annotate DMRs
# TODO: Load dependencies outside the function, specially those shared with other functions?
# TODO: Add loop in both replicate and noReplicate that go over every pair and produce tables and plots for each one.
# TODO: Which objects we want the function to return?
annotateDMRs<- function(dmrs.gr,chrsOfInterest=paste("chr",1:22,sep=""),genome.version='hg19') {
  # Dependencies
  suppressMessages(library('annotatr'))
  suppressMessages(library('Biobase'))
  suppressMessages(library('AnnotationDbi'))
  suppressMessages(library('GenomicFeatures'))
  suppressMessages(library('IRanges'))
  suppressMessages(library('rtracklayer'))
  suppressMessages(library('GenomicRanges'))
  suppressMessages(library('ggplot2'))
  suppressMessages(library('TxDb.Hsapiens.UCSC.hg19.knownGene'))
  suppressMessages(library('Homo.sapiens'))

  # Complete annotations
  write(paste("[",date(),"]: Building annotations..."), stdout())
  annots <- c('_genes_promoters','_genes_1to5kb','_genes_cds','_genes_5UTRs','_genes_exons','_genes_firstexons','_genes_introns',
	      '_genes_intronexonboundaries','_genes_exonintronboundaries','_genes_3UTRs','_genes_intergenic','_cpg_islands',
	      '_cpg_shores','_cpg_shelves','_cpg_inter','_enhancers_fantom','_basicgenes','_cpgs')
  annots <- paste(genome.version,annots,sep="")
  annots.all.gr <- build_annotations(genome=genome.version,annotations=annots)
  
  # Subset gene body for abbreviated version and get rid of those annotations with SYMBOL=<NA>
  annots <- c('_genes_promoters','_genes_cds','_genes_5UTRs','_genes_exons','_genes_firstexons','_genes_introns',
	      '_genes_intronexonboundaries','_genes_exonintronboundaries','_genes_3UTRs')
  annots <- paste(genome.version,annots,sep="")
  annots.genes.gr <- annots.all.gr[annots.all.gr$type %in% annots]
  annots.genes.gr <- annots.all.gr[!is.na(annots.all.gr$symbol),]

  # Assign ID to each DMR
  dmrs.gr$dmrId <- seq(1,length(dmrs.gr))

  # Annotate DMRs with annotatr for plots with annotatr
  write(paste("[",date(),"]: Default annotation with annotatr..."), stdout())
  dmrs.annotatr.gr <- annotate_regions(dmrs.gr,annots.all.gr,ignore.strand=TRUE)

  # Annotate DMR (overlapping) with detail manually
  # Note: this annotation agrees with annotatr's but contains metadata that can be included in the output table
  # Note: all the annotations are included here (gene+enhancers+etc)
  write(paste("[",date(),"]: Detailed annotation for overlapping DMRs..."), stdout())
  olap.gr <- findOverlaps(dmrs.gr,annots.all.gr,ignore.strand=TRUE)
  id.list <- strsplit(as.character(annots.all.gr[subjectHits(olap.gr),]$id),":")
  id.df <- do.call(rbind,id.list)
  dmrs.olap.gr <- dmrs.gr[queryHits(olap.gr),]
  dmrs.olap.gr$feature <- id.df[,1]
  dmrs.olap.gr$feature_id <- id.df[,2]
  dmrs.olap.gr$tx_id <- annots.all.gr[subjectHits(olap.gr),]$tx_id
  dmrs.olap.gr$symbol <- annots.all.gr[subjectHits(olap.gr),]$symbol
  dmrs.olap.df <- as.data.frame(dmrs.olap.gr)
  dmrs.olap.df$strand <- NULL
  dmrs.olap.df$width <- NULL
  
  # Sort DMRs based on # features detected
  dmrs.olap.table <- as.data.frame(table(dmrs.olap.gr$dmrId))
  colnames(dmrs.olap.table) <- c('dmrId','Freq')
  dmrs.olap.table <- dmrs.olap.table[order(dmrs.olap.table$Freq,decreasing=TRUE),]
  dmrs.olap.df$rkg <- match(dmrs.olap.df$dmrId,dmrs.olap.table$dmrId)
  dmrs.olap.df <- dmrs.olap.df[order(dmrs.olap.df$rkg),]

  # Annotate DMR based only on gene body and collapse
  write(paste("[",date(),"]: Abbreviated gene body annotation..."), stdout())
  
  # (1) Overlapping
  # Note: for precede and follows the DMRs found here are excluded
  olap.genes.gr <- findOverlaps(dmrs.gr,annots.genes.gr,select="all",ignore.strand=TRUE)
  dmrs.olap.genes.gr <- dmrs.gr[queryHits(olap.genes.gr),]
  dmrs.olap.genes.gr$symbol <- annots.genes.gr[subjectHits(olap.genes.gr)]$symbol 
  dmrs.olap.genes.gr <- unique(as.data.frame(dmrs.olap.genes.gr))
  dmrs.olap.genes.df <- as.data.frame(dmrs.olap.genes.gr)
  dmrs.olap.genes.df$distance <- rep(0,nrow(dmrs.olap.genes.df))
  dmrs.olap.genes.df$strand <- NULL
  dmrs.olap.genes.df$width <- NULL
  dmrs.olap.genes.df <- dmrs.olap.genes.df[!is.na(dmrs.olap.genes.df$symbol),]
  dmrs.olap.found.df <- unique(sort(dmrs.olap.genes.df$dmrId))
  
  # (2) DMR precedes gene (left->right)
  # Note: here overlapping DMRs are NOT excluded and genes overlapping can appear here again because the overlap between DMR
  # and gene is not total
  # Note: distance has positive sign since the gene follows the DMR
  prec.genes.gr <- precede(dmrs.gr,annots.genes.gr,select="all",ignore.strand=TRUE)
  dmrs.prec.genes.gr <- dmrs.gr[queryHits(prec.genes.gr),]
  dmrs.prec.genes.gr$symbol <- annots.genes.gr[subjectHits(prec.genes.gr)]$symbol
  dmrs.prec.genes.gr$distance <- distance(dmrs.prec.genes.gr,annots.genes.gr[subjectHits(prec.genes.gr)])
  dmrs.prec.genes.gr <- unique(as.data.frame(dmrs.prec.genes.gr))
  dmrs.prec.genes.df <- as.data.frame(dmrs.prec.genes.gr)
  dmrs.prec.genes.df$strand <- NULL
  dmrs.prec.genes.df$width <- NULL
  dmrs.prec.genes.df <- dmrs.prec.genes.df[!is.na(dmrs.prec.genes.df$symbol),]
  # dmrs.prec.genes.df <- dmrs.prec.genes.df[!(dmrs.prec.genes.df$dmrId %in% dmrs.olap.found.df),]
  
  # (3) DMR follows gene (left->right)
  # Note: here overlapping DMRs are NOT excluded and genes overlapping can appear here again because the overlap between DMR
  # and gene is not total
  # Note: distance has negative sign since the gene precedes the DMR
  foll.genes.gr <- follow(dmrs.gr,annots.genes.gr,select="all",ignore.strand=TRUE)
  dmrs.foll.genes.gr <- dmrs.gr[queryHits(foll.genes.gr),]
  dmrs.foll.genes.gr$symbol <- annots.genes.gr[subjectHits(foll.genes.gr)]$symbol 
  dmrs.foll.genes.gr$distance <- -distance(dmrs.foll.genes.gr,annots.genes.gr[subjectHits(foll.genes.gr)])
  dmrs.foll.genes.gr <- unique(as.data.frame(dmrs.foll.genes.gr))
  dmrs.foll.genes.df <- as.data.frame(dmrs.foll.genes.gr)
  dmrs.foll.genes.df$strand <- NULL
  dmrs.foll.genes.df$width <- NULL
  dmrs.foll.genes.df <- dmrs.foll.genes.df[!is.na(dmrs.foll.genes.df$symbol),]
  # dmrs.foll.genes.df <- dmrs.foll.genes.df[!(dmrs.foll.genes.df$dmrId %in% dmrs.olap.found.df),]
  
  # Collapse the searches and sort DMRs based on # features detected in complete annotation (beyond gene body)
  dmrs.all.genes.df <- rbind(dmrs.olap.genes.df,dmrs.prec.genes.df,dmrs.foll.genes.df)
  dmrs.all.genes.df$rkg <- match(dmrs.all.genes.df$dmrId,dmrs.olap.table$dmrId)
  dmrs.all.genes.df <- dmrs.all.genes.df[order(dmrs.all.genes.df$rkg,dmrs.all.genes.df$distance),]

  # Annotate DMR based only on TSS and collapse
  write(paste("[",date(),"]: Abbreviated TSS annotation..."), stdout())
  
  # Get TSS from TxDb
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  chr.filter <- list(tx_chrom = chrsOfInterest)
  annots.tss.gr <- promoters(transcripts(txdb,filter=chr.filter),upstream=0L,downstream=1L)
  symbol.table <- select(Homo.sapiens,key=as.character(annots.tss.gr$tx_id),keytype="TXID",columns=c("SYMBOL"))
  annots.tss.gr$symbol <- symbol.table$SYMBOL
  annots.tss.gr <- unique(annots.tss.gr[!is.na(annots.tss.gr$symbol),])
  
  # (1) Overlapping
  # Note: for precede and follows the DMRs found here are excluded
  olap.tss.gr <- findOverlaps(dmrs.gr,annots.tss.gr,select="all",ignore.strand=TRUE)
  dmrs.olap.tss.gr <- dmrs.gr[queryHits(olap.tss.gr),]
  dmrs.olap.tss.gr$symbol <- annots.tss.gr[subjectHits(olap.tss.gr)]$symbol 
  dmrs.olap.tss.gr <- unique(as.data.frame(dmrs.olap.tss.gr))
  dmrs.olap.tss.df <- as.data.frame(dmrs.olap.tss.gr)
  dmrs.olap.tss.df$distance <- rep(0,nrow(dmrs.olap.tss.df))
  dmrs.olap.tss.df$strand <- NULL
  dmrs.olap.tss.df$width <- NULL
  dmrs.olap.tss.df <- dmrs.olap.tss.df[!is.na(dmrs.olap.tss.df$symbol),]
  
  # (2) DMR precedes TSS (left->right)
  # Note: distance has positive sign since the gene follows the DMR
  prec.tss.gr <- precede(dmrs.gr,annots.tss.gr,select="all",ignore.strand=TRUE)
  dmrs.prec.tss.gr <- dmrs.gr[queryHits(prec.tss.gr),]
  dmrs.prec.tss.gr$symbol <- annots.tss.gr[subjectHits(prec.tss.gr)]$symbol
  dmrs.prec.tss.gr$distance <- distance(dmrs.prec.tss.gr,annots.tss.gr[subjectHits(prec.tss.gr)])
  dmrs.prec.tss.gr <- unique(as.data.frame(dmrs.prec.tss.gr))
  dmrs.prec.tss.df <- as.data.frame(dmrs.prec.tss.gr)
  dmrs.prec.tss.df$strand <- NULL
  dmrs.prec.tss.df$width <- NULL
  dmrs.prec.tss.df <- dmrs.prec.tss.df[!is.na(dmrs.prec.tss.df$symbol),]
  
  # (3) DMR follows TSS (left->right)
  # Note: distance has negative sign since the gene precedes the DMR
  foll.tss.gr <- follow(dmrs.gr,annots.tss.gr,select="all",ignore.strand=TRUE)
  dmrs.foll.tss.gr <- dmrs.gr[queryHits(foll.tss.gr),]
  dmrs.foll.tss.gr$symbol <- annots.tss.gr[subjectHits(foll.tss.gr)]$symbol 
  dmrs.foll.tss.gr$distance <- -distance(dmrs.foll.tss.gr,annots.tss.gr[subjectHits(foll.tss.gr)])
  dmrs.foll.tss.gr <- unique(as.data.frame(dmrs.foll.tss.gr))
  dmrs.foll.tss.df <- as.data.frame(dmrs.foll.tss.gr)
  dmrs.foll.tss.df$strand <- NULL
  dmrs.foll.tss.df$width <- NULL
  dmrs.foll.tss.df <- dmrs.foll.tss.df[!is.na(dmrs.foll.tss.df$symbol),]
  
  # Collapse the searches and sort DMRs based on # features detected in complete annotation (beyond gene body)
  dmrs.all.tss.df <- rbind(dmrs.olap.tss.df,dmrs.prec.tss.df,dmrs.foll.tss.df)
  dmrs.all.tss.df$rkg <- match(dmrs.all.tss.df$dmrId,dmrs.olap.table$dmrId)
  dmrs.all.tss.df <- dmrs.all.tss.df[order(dmrs.all.tss.df$rkg,dmrs.all.tss.df$distance),]

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
