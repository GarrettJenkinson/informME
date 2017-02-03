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

# FastaToCpG.sh

# last modified 12/09/16

# Shell script to run function FastaToCpG.m in MATLAB

# Get inputs
FastaFile=$1
species=$2
MATLICE=$3

echo "starting command: FastaToCpG('${FastaFile}','species','${species}');"

# Run MATLAB command
matlab -nodesktop -singleCompThread -nosplash -nodisplay -c ${MATLICE} -r "disp('Job Running');tic;FastaToCpG('${FastaFile}','species','${species}');toc;"

# Check if error in MATLAB, otherwise declare success
EXITCODE=$?
if [ $EXITCODE -ne 0 ]
then
   echo "error in MATLAB"
   echo "tic;FastaToCpG('${FastaFile}','species','${species}');toc;"
   exit $EXITCODE
fi

echo "command successful:"
echo "tic;FastaToCpG('${FastaFile}','species','${species}');toc;"
