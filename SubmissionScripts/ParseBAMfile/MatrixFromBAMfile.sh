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

# MatrixFromBAMfile.sh

# last modified 12/09/16

# Shell script to run function MatrixFromBAMfile.m in MATLAB

# Get inputs
bamFilename=$1
chr_num=$2
totalProcessors=$3 
processorNum=$4
species=$5
includeChrInRef=$6
numBasesToTrim=$7
MATLICE=$8
BAMPATH=$9

echo "starting command: MatrixFromBAMfile('${bamFilename}', ${chr_num}, 'totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','includeChrInRef',${includeChrInRef},'numBasesToTrim',${numBasesToTrim},'bamFilePathRoot','${BAMPATH}');"

# Run MATLAB command
matlab -singleCompThread -nosplash -nodisplay -c ${MATLICE} -r "disp('Job Running');tic; MatrixFromBAMfile('${bamFilename}', ${chr_num}, 'totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','includeChrInRef',${includeChrInRef},'numBasesToTrim',${numBasesToTrim},'bamFilePathRoot','${BAMPATH}');toc"

# Check if error in MATLAB, otherwise declare success
EXITCODE=$?
if [ $EXITCODE -ne 0 ]
then
   echo "error in MATLAB"
   echo "tic; 
   exit $EXITCODE
fi

echo "command successful"
echo "tic; MatrixFromBAMfile('${bamFilename}', ${chr_num}, 'totalProcessors',${totalProcessors},'processorNum',${processorNum},'species','${species}','includeChrInRef',${includeChrInRef},'numBasesToTrim',${numBasesToTrim},'bamFilePathRoot','${BAMPATH}');toc"
