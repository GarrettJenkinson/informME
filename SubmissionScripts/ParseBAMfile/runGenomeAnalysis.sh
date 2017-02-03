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

# runGenomeAnalysis.sh

# last modified 12/09/16

# Main script to run the script FastaToCpG.sh on the 
# Sun Grid Engine cluster

# Get inputs
FastaFile=$1
species=$2

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"

# Non-user specified directories
# (trailing slash)
SCRIPTDIR="$IPATH/SubmissionScripts/ParseBAMfile/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/ParseBAMfile/"

# Run FastaToCpG.sh script
JOB=`qsub << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_GA 
#$ -l mem_free=15G
#$ -l h_vmem=16G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}FastaToCpG.sh "${FastaFile}" "${species}" "${MATLICE}"
EOJ
`
echo "${JOB} for file ${INFILENAME3} submitted on `date`"

