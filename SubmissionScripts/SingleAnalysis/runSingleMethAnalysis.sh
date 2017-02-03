#!/bin/bash
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

# runSingleMethAnalysis.sh

# last modified 12/09/16

# Main script to run the scripts SingleMethAnalysis.sh, 
# MergeSingleMethAnalysis.sh, and MakeSingleMethAnalysisBEDs.sh 
# on the Sun Grid Engine cluster

# Get inputs
phenoName=$1
species=$2

# USER SPECIFIED VARIABLES
ESIflag=0 # Set flag to 1 if entropy sensitivity computation is desired.
MCflag=0  # Set flag to 1 if methylation channel computations are desired.

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"

# Non-user specified directories
# (trailing slash)
SCRIPTDIR="$IPATH/SubmissionScripts/SingleAnalysis/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/SingleAnalysis/"

# Specify total processors per chromosome
totalProcessors=50

# Set values based on input species
if [ "${species}" = "Chicken" ]; then
    totalChrNum=28
elif [ "${species}" = "Mouse" ]; then
    totalChrNum=19
elif [ "${species}" = "Opossum" ]; then
    totalChrNum=8
else # assume its "Human"
    totalChrNum=22
fi

JOBLIST2=""

# Perform single methylation analysis on each chromosome
for ((chr=1; chr<=totalChrNum; chr++))
do # start do over chromosomes
JOBLIST1=""
for ((processorNum=1; processorNum<=totalProcessors; processorNum++))
do # start do over processors
    
# Run SingleMethAnalysis.sh script
JOB=`qsub << EOJ
#
#$ -wd $WORKDIR
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_SMA
#$ -l mem_free=20G
#$ -l h_vmem=21G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}SingleMethAnalysis.sh "$chr" "$phenoName" "$totalProcessors" "$processorNum" "$species" "$ESIflag" "$MCflag" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"
if [ "${JOBLIST1}" = "" ]; then
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
    JOBLIST1="${TEMP}"
else
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
    JOBLIST1="${JOBLIST1},${TEMP}"
fi
done # end of loop over processors

# # Run MergeSingleMethAnalysis.sh script
JOB=`qsub -hold_jid ${JOBLIST1} << EOJ
#
#$ -wd $WORKDIR
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_SMA_MERGE
#$ -l mem_free=61G
#$ -l h_vmem=62G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
#
${SCRIPTDIR}MergeSingleMethAnalysis.sh "$chr" "$phenoName" "$species" "$totalProcessors" "$ESIflag" "$MCflag" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"
if [ "${JOBLIST2}" = "" ]; then
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
    JOBLIST2="${TEMP}"
else
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
    JOBLIST2="${JOBLIST2},${TEMP}"
fi
done # end of do over chromosomes

# Run MakeSingleMethAnalysisBEDs.sh script
JOB=`qsub -hold_jid ${JOBLIST2} << EOJ
#
#$ -wd $WORKDIR
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_SMA_MBED
#$ -l mem_free=61G
#$ -l h_vmem=62G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
#
${SCRIPTDIR}MakeSingleMethAnalysisBEDs.sh "$phenoName" "$species" "1" "$totalChrNum" "$ESIflag" "$MCflag" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"
