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

# runDiffMethAnalysis.sh

# last modified 12/09/16

# Main script to run script MakeDiffMethAnalysisBEDs.sh on the 
# Sun Grid Engine cluster. 

# Get inputs
phenoName_t=$1
phenoName_r=$2
species=$3

# USER SPECIFIED VARIABLES
ESIflag=0 # Set flag to 1 if entropy sensitivity computation is desired.
MCflag=0  # Set flag to 1 if methylation channel computations are desired.

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"

# Non-user specified directories
# (trailing slash)
SCRIPTDIR="$IPATH/SubmissionScripts/DiffAnalysis/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/DiffAnalysis/"

# Set values based on input species
if [ "${species}" = "Chicken" ]; then
    max_chr_num=28
elif [ "${species}" = "Mouse" ]; then
    max_chr_num=19
elif [ "${species}" = "Opossum" ]; then
    max_chr_num=8
else # assume its "Human"
    max_chr_num=22
fi

# Run MakeDiffMethAnalysisBEDs.sh script
JOB=`qsub << EOJ
#
#$ -wd $WORKDIR
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_DMA
#$ -l mem_free=21G
#$ -l h_vmem=22G
#$ -l h_fsize=100G
#$ -l cegs2
#$ -q cegs2.q
#$ -o ${SOUTDIR}
#
${SCRIPTDIR}MakeDiffMethAnalysisBEDs.sh  "{'$phenoName_t','$phenoName_r'}" "$species" "1" "$max_chr_num" "$ESIflag" "$MCflag" "$MATLICE"
EOJ
`
echo "JobID = ${JOB} submitted on `date`"
TEMP=$(echo "${JOB}" | grep -Eo '[0-9]{1,7}')
JOBLIST="${TEMP}"
