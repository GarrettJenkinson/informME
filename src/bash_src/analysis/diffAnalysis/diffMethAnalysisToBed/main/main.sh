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

# Reference from fastaToCpg
prefix_1="toy_normal"
prefix_2="toy_cancer"
estimation="../../../../modeling/mergeEstimation/main/out/"
analysis="../../../singleAnalysis/mergeSingleMethAnalysis/main/out/"
reference="../../../../parseBamFile/fastaToCpg/main/"

# Execute command
../diffMethAnalysisToBed.sh -r "$reference" --analdir_1 "$analysis" --analdir_2 "$analysis" -d out --min_chr 1 --max_chr 5 --MC --ESI -- "$prefix_1" "$prefix_2" 1>/dev/null
