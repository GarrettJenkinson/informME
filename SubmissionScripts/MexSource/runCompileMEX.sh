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

# runCompileMEX.sh

# last modified 12/09/16

# Main script to run the script compile.sh on the 
# Sun Grid Engine cluster

# USER SPECIFIED DIRECTORIES
# (no trailing slash)
IPATH="/path/to/informME_folder"
MATLICE="/path/to/MATLAB_license_directory/[filename].lic"
PATH1="/path/to/eigen_instalation"
PATH2="/path/to/mpfr_instalation"

# Non-user specified directories
# (trailing slash)
EIGEN="$PATH1/"
MPFRINCL="$PATH2/include/"
MPFRLIB="$PATH2/lib/"
SCRIPTDIR="$IPATH/SubmissionScripts/MexSource/"
SOUTDIR="$IPATH/SubmissionScripts/STDouts/"
WORKDIR="$IPATH/MexSource/"

# Run compile.sh script
JOB=`qsub << EOJ
#
#$ -wd ${WORKDIR}
#$ -j y
#$ -r y
#$ -S /bin/bash
#$ -N INFM_CMP
#$ -l mem_free=15G
#$ -l h_vmem=16G
#$ -l h_fsize=100G
#$ -o ${SOUTDIR}
# 
${SCRIPTDIR}compile.sh "${EIGEN}" "${MPFRINCL}" "${MPFRLIB}" "${MATLICE}"
EOJ
`
echo "${JOB} for COMPILE submitted on `date`"
