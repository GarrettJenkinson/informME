#!/bin/bash
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# or see <http://www.gnu.org/licenses/>.

# runModelEstimation.sh

# last modified 12/09/16

# Main script to run the scripts Estimation.sh and 
# MergeEstimation.sh on the Sun Grid Engine cluster

# Get inputs
bamFileNames=$1
phenoName=$2
species=$3

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"

# Non-user specified directories
# (trailing slash)
SCRIPTDIR="$IPATH/SubmissionScripts/Modeling/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/Modeling/"

# Specify total processors per chromosome
totalProcessors=200

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

# Process BAM files to generate estimation results
for ((chr=1; chr<=totalChrNum; chr++))
do # start do over chromosomes
JOBLIST=""
for ((processorNum=1; processorNum<=totalProcessors; processorNum++))
do # start do over processors
   JOB=`qsub << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_EST 
#$ -l mem_free=12G
#$ -l h_vmem=13G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}Estimation.sh "$bamFileNames" "$chr" "$phenoName" "$species" "$totalProcessors" "$processorNum" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"
if [ "${JOBLIST}" = "" ]; then
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
	JOBLIST="${TEMP}"
else
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
	JOBLIST="${JOBLIST},${TEMP}"
fi
done # end of loop over processors

# Merge estimation results
JOB=`qsub -hold_jid ${JOBLIST} << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_EST_MERGE
#$ -l mem_free=61G
#$ -l h_vmem=62G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}MergeEstimation.sh "$bamFileNames" "$chr" "$phenoName" "$species" "$totalProcessors" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"

done # end of do over chromosomes
