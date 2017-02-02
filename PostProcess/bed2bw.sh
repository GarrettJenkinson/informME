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



# bed2bw

# last modified 12/09/16

# Script to convert BED files to BigWig files

# Get inputs

BEDDIR=$1
ASSEMBLY=$2

# Make list of BED files in BEDDIR
FILES="${BEDDIR}*.bed"

cd ${BEDDIR}
mkdir BWfiles

# Read chr sizes
fetchChromSizes ${ASSEMBLY} > ${ASSEMBLY}.sizes

# Loop over files
for f in $FILES
do
  echo "Processing $f file..."
  tail -n +2 $f | awk -v OFS='\t' '{print $1,$2,$3,$4}' > $f.nh
  bedtools slop -i $f.nh -g ${ASSEMBLY}.sizes -b 0 | bedClip stdin ${ASSEMBLY}.sizes $f.cl
  LC_COLLATE=C sort -k1,1 -k2,2n $f.cl > $f.st
  bedGraphToBigWig $f.st ${ASSEMBLY}.sizes ${f/bed/bw}
  rm $f.nh
  rm $f.cl
  rm $f.st
done

# Clean up:
rm ${BEDDIR}${ASSEMBLY}.sizes
mv *.bw BWfiles/

