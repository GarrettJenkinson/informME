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

# Extend regex capabilities
shopt -s extglob

# Get directories
abspath_script="$(readlink -f -e "$0")"
script_absdir="$(dirname "$abspath_script")"
script_name="$(basename "$0" .sh)"

# Print help if no arguments
if [ $# -eq 0 ]
    then
        cat "$script_absdir/${script_name}_help.txt"
        exit 1
fi

# Getopt command
TEMP=$(getopt -q -o hd: -l help,outdir: -n "$script_name.sh" -- "$@")

if [ $? -ne 0 ]
then
  echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
  cat "$script_absdir"/${script_name}_help.txt
  echo "[$(date)]: Terminating" >&2
  exit -1
fi

eval set -- "$TEMP"

# Defaults
outdir="$PWD"

# Options
while true
do
  case "$1" in
    -h|--help)
      cat "$script_absdir"/${script_name}_help.txt
      exit
      ;;  
    -d|--outdir)
      outdir="$2"
      shift 2
      ;;  
    --) 
      shift
      break
      ;;  
    *)
      echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
      cat "$script_absdir"/${script_name}_help.txt
      echo "[$(date)]: Terminating" >&2
      exit -1
      ;;
  esac
done

# Check number of arguments and copy them
if [ "$#" -ne 2 ]; then
   echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
   cat "$script_absdir"/${script_name}_help.txt
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
bed_dir="$1"
bed_path="$(readlink -f -e "$bed_dir")"
genome="$2"

# Check valid outdir
if [ -z "$outdir" ];then
   echo -e "[$(date)]: \e[31mERROR: Output directory is empty string ...\e[0m" >&2
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
mkdir -p "$outdir"

# Start time
SECONDS=0
echo "[$(date)]: Starting ..." 

# Fetch chromosome sizes
echo "[$(date)]: Fetching chromosome sizes ..." 
temp_chr_sizes="${outdir}/${genome}.chrom.sizes"
fetchChromSizes $genome > "$temp_chr_sizes"

# Loop over BED files in bed_dir
for f in "$bed_path"/*.bed
do
	echo "[$(date)]: Processing ${f} ..." 
    filename="$(basename "$f" .bed)"
    outfile="${outdir}/${filename}.bw"
  	tail -n +2 $f | awk -v OFS='\t' '{print $1,$2,$3,$4}' > "${f}.nohead"
  	bedtools slop -i "$f.nohead" -g "$temp_chr_sizes" -b 0 | bedClip stdin "$temp_chr_sizes" "${f}.clip"
  	LC_COLLATE=C sort -k1,1 -k2,2n "${f}.clip" > "${f}.sorted"
  	bedGraphToBigWig "${f}.sorted" "$temp_chr_sizes"  "$outfile"
  	rm "${f}.nohead"
  	rm "${f}.clip"
  	rm "${f}.sorted"
done

# Clean
rm "$temp_chr_sizes"

# Time elapsed
time_elapsed="$SECONDS"
echo "[$(date)]: Total time elapsed: $(( $time_elapsed / 3600)) h \
$(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
