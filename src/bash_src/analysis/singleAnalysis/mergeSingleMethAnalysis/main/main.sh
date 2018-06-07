#!/usr/bin/env bash
#
# informME: An information-theoretic pipeline for WGBS data
# Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu),
# and Jordi Abante (jabante1@jhu.edu)
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
#

# Remove previous results
rm -rf out
cp -r ../../singleMethAnalysis/main/out/ input

# Execute command for normal toy
echo "Processing normal..."
seq 5 | xargs -I {X} \
bash -c '../mergeSingleMethAnalysis.sh -r ../../../../parseBamFile/fastaToCpg/main/ -e ../../../../modeling/mergeEstimation/main/out/ -a input/ -d out/ --MC --ESI --MSI toy_normal {X} 4' 1>/dev/null

# Execute command for cancer toy
echo "Processing cancer..."
seq 5 | xargs -I {X} \
bash -c '../mergeSingleMethAnalysis.sh -r ../../../../parseBamFile/fastaToCpg/main/ -e ../../../../modeling/mergeEstimation/main/out/ -a input/ -d out/ --MC --ESI --MSI toy_cancer {X} 4' 1>/dev/null

# Execute command for pooled toy
echo "Processing pooled..."
seq 5 | xargs -I {X} \
bash -c '../mergeSingleMethAnalysis.sh -r ../../../../parseBamFile/fastaToCpg/main/ -e ../../../../modeling/mergeEstimation/main/out/ -a input/ -d out/ --MC --ESI --MSI toy_pooled {X} 4' 1>/dev/null

# Remove input
rm -rf input
