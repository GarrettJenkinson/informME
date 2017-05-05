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

# runDataMatrixGeneration.sh

# last modified 12/09/16

# Main script to run the scripts MatrixFromBAMfile.sh 
# and MergeMatrices.sh on the Sun Grid Engine cluster

# Get inputs
FILENAME=$1
species=$2
numBasesToTrim=$3

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"

# Non-user specified directories
# (trailing slash)
BAMPATH="$IPATH/ParseBAMfile/indexedBAMfiles/"
SCRIPTDIR="$IPATH/SubmissionScripts/ParseBAMfile/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/ParseBAMfile/"

# Specify total processors per chromosome
totalProcessors=200

# Set values based on input species
if [ "${species}" = "Chicken" ]; then
    includeChr=1
    totChrNum=28
elif [ "${species}" = "Mouse" ]; then
    includeChr=1
    totChrNum=19
elif [ "${species}" = "Opossum" ]; then
    includeChr=1
    totChrNum=8
else # assume its Human
    includeChr=0
    totChrNum=22
fi

# Process BAM file to generate methylation matrices
for ((chr=1; chr<=totChrNum; chr++))
do # start do over chromosomes
JOBLIST=""
for ((processorNum=1; processorNum<=totalProcessors; processorNum++))
do # start do over processors

# Run MatrixFromBAMfile.sh script
JOB=`qsub << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_MTRX 
#$ -l mem_free=6G
#$ -l h_vmem=7G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}MatrixFromBAMfile.sh "$FILENAME" "$chr" "$totalProcessors" "$processorNum" "$species" "$includeChr" "$numBasesToTrim" "$MATLICE" "$BAMPATH"
EOJ
`
echo "${JOB} for file ${FILENAME} submitted on `date`"
if [ "${JOBLIST}" = "" ]; then
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,9}')
	JOBLIST="${TEMP}"
else
    TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,9}')
	JOBLIST="${JOBLIST},${TEMP}"
fi
done # end do over processors

# Run MergeMatrices.sh script
JOB=`qsub -hold_jid ${JOBLIST} << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_MTRX_MERGE
#$ -l mem_free=31G
#$ -l h_vmem=32G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}MergeMatrices.sh "$FILENAME" "$chr" "$totalProcessors" "$species" "$includeChr" "$numBasesToTrim" "$MATLICE" "$BAMPATH"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"

done # end do over chromosomes
