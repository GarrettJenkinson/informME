#!/bin/bash
#
# informME: An information-theoretic pipeline for WGBS data
# Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# or see <http://www.gnu.org/licenses/>.

# SingleMethAnalysis.sh

# last modified 12/09/16

# Shell acript to run function MethAnalysisForChr.m in MATLAB

# get inputs
chr_num=$1
phenoName=$2
totalProcessors=$3
processorNum=$4
species=$5
ESIflag=$6
MCflag=$7
MATLICE=$8

echo "starting command: MethAnalysisForChr(${chr_num},'${phenoName}','totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','ESIflag',${ESIflag},'MCflag',${MCflag});"

# Run MATLAB command

matlab -nodesktop -singleCompThread -nosplash -nodisplay -c ${MATLICE} -r "try, disp('Job Running');tic;MethAnalysisForChr(${chr_num},'${phenoName}','totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','ESIflag',${ESIflag},'MCflag',${MCflag});toc, catch, exit(1), end, exit(0);"

# Check if error in MATLAB, otherwise declare success
EXITCODE=$?
if [ $EXITCODE -ne 0 ]
then
   echo "error in MATLAB"
   echo "tic;MethAnalysisForChr(${chr_num},'${phenoName}','totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','ESIflag',${ESIflag},'MCflag',${MCflag});toc;"
   exit $EXITCODE
fi

echo "command successful"
echo "tic;MethAnalysisForChr(${chr_num},'${phenoName}','totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','ESIflag',${ESIflag},'MCflag',${MCflag});toc;"
