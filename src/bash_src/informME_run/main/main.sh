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
rm -rf intermediate
rm -rf tmp

# Process normal
echo "Processing normal..."
seq 5 | xargs -I {X} --max-proc 1 bash -c '../informME_run.sh -r ../../parseBamFile/fastaToCpg/main/ -e intermediate --tmp tmp -d out -q 4 --time_limit 2 --total_part 4 -- toy_normal_pe toy_normal {X} 4' 1>/dev/null

# Process cancer
echo "Processing cancer..."
seq 5 | xargs -I {X} --max-proc 1 bash -c '../informME_run.sh -r ../../parseBamFile/fastaToCpg/main/ -e intermediate --tmp tmp -d out -q 4 --time_limit 2 --total_part 4 -- toy_cancer_pe toy_cancer {X} 4' 1>/dev/null

# Process pooled
echo "Processing pooled..."
seq 5 | xargs -I {X} --max-proc 1 bash -c '../informME_run.sh -r ../../parseBamFile/fastaToCpg/main/ -e intermediate --tmp tmp -d out -q 4 --time_limit 2 --total_part 4 -- toy_normal_pe,toy_cancer_pe toy_pooled {X} 4' 1>/dev/null
