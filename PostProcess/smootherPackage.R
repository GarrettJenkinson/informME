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


# The following libraries are required:
# library(rtracklayer)
# library(logitnorm)
# library(mixtools)

# THIS FILE WILL PROVIDE EXAMPLE USAGE OF THE jsDMR.R file's functions

rm(list=ls())

#
#######################################################################
# Initializations
#######################################################################
#

# define input and output directories
inFolder <- "/home/gar/Desktop/BWfiles/"
outFolder <- "/home/gar/Desktop/Rscripts/Smoother/Bioinformatics/"

# null comparisons
file1 <- "JSD-lungnormal12-VS-lungnormal10.bed"
file2 <- "JSD-lungnormal12-VS-lungnormal8.bed"
file3 <- "JSD-lungnormal8-VS-lungnormal10.bed"
nullFiles <- c(file1,file2,file3)
  
# test cases
file4 <- "JSD-lungcancer7-VS-lungnormal8.bed"
file5 <- "JSD-lungcancer9-VS-lungnormal10.bed"
file6 <- "JSD-lungcancer11-VS-lungnormal12.bed"
altFiles <- c(file4,file5,file6)

# pooled case. might not work correctly with null model, 
# but should work in non-replicate mixture case
file7 <- "JSD-lungcancerP-VS-lungnormalP.bed"

#
#######################################################################
# Function Definitions
#######################################################################
#

source(paste(outFolder,"jsDMR.R",sep=""))

#
#######################################################################
# Do empirical null modeling
#######################################################################
#

DMRs <- runEmpiricalNullDMR(nullFiles,altFiles,inFolder,outFolder)

#
#######################################################################
#  Do gaussian mixture modeling for non-replicate case
#######################################################################
#

DMRmix <- runMixLogitNullDMR(file4,inFolder,outFolder)
DMRmix <- runMixLogitNullDMR(file5,inFolder,outFolder)
DMRmix <- runMixLogitNullDMR(file6,inFolder,outFolder)
DMRmix <- runMixLogitNullDMR(file7,inFolder,outFolder)
