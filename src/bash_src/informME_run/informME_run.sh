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

# Get default locations
INTERDIR="$PWD"
if ! [ -x "$(command -v informME_source_config)" ]; then
   # If bin folder not in path, exit with Error
   echo -e "[$(date)]: \e[31mERROR: bin folder not in PATH. Exiting with error ...\e[0m" >&2
   exit 1
else
   # If bin folder added to path, source config file
   source informME_source_config 1>/dev/null
fi

# Getopt command
TEMP="$(getopt -q -o hr:m:e:d:l:q: -l help,refdir:,matdir:,estdir:,tmpdir:,outdir:,MATLICENSE:,threads:,time_limit:,total_part: -n "$script_name.sh" -- "$@")"

if [ $? -ne 0 ]
then
  echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
  echo "[$(date)]: Terminating" >&2
  exit -1
fi

eval set -- "$TEMP"

# Defaults
refdir="$REFGENEDIR"
matdir="$INTERDIR"
estdir="$INTERDIR"
tmpdir="$SCRATCHDIR"
outdir="$INTERDIR"
threads=1
time_limit=60
total_part=200

# Options
while true
do
  case "$1" in
    -h|--help)
      cat "$script_absdir"/${script_name}_help.txt
      exit
      ;;  
    -r|--refdir)
      refdir="$2"
      shift 2
      ;;
    -m|--matdir)
      matdir="$2"
      shift 2
      ;;
    -e|--estdir)
      estdir="$2"
      shift 2
      ;;
    --tmpdir)
      tmpdir="$2"
      shift 2
      ;;
    -d|--outdir)
      outdir="$2"
      shift 2
      ;;
    -l|--MATLICENSE)
      MATLICENSE="$2"
      shift 2
      ;;  
    -q|--threads)
      threads="$2"
      shift 2
      ;;  
    --time_limit)
      time_limit="$2"
      shift 2
      ;;
    --total_part)
      total_part="$2"
      shift 2
      ;;
    --) 
      shift
      break
      ;;  
    *)
      echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
      echo "[$(date)]: Terminating" >&2
      exit -1
      ;;
  esac
done

# Check number of arguments and copy them
if [ "$#" -ne 3 ]; then
   echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
mat_files="$1"
prefix="$2"
chr_num="$3"

# Check valid outdir
if [ -z "$outdir" ];then
   echo -e "[$(date)]: \e[31mERROR: Output directory is empty string ...\e[0m" >&2
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
mkdir -p "$outdir"
mkdir -p "${outdir}/chr${chr_num}"

#check if output already exists
if [ -r "${outdir}/chr${chr_num}/${prefix}_analysis.mat" ]
then
  echo -e "[$(date)]:\e[31m Warning: Final output file:\e[0m" >&2
  echo -e "[$(date)]:\e[31m ${outdir}/chr${chr_num}/${prefix}_analysis.mat\e[0m" >&2
  echo -e "[$(date)]:\e[31m already exists. Delete file and rerun if you wish to overwrite.\e[0m" >&2 
  echo -e "[$(date)]:\e[31m Terminating...\e[0m" >&2
  exit 0
fi

# Run estimation on chunks
SECONDS=0
echo "[$(date)]: Call: estimation.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
seq "$total_part" | xargs -I {X} --max-proc "$threads" bash -c "timeout --signal=SIGINT '$time_limit'm estimation.sh -r '$refdir' -m '$matdir' -d '$tmpdir' -- '$mat_files' '$prefix' '$chr_num' '$total_part' {X}"

# Check if everything OK
EXITCODE="$?"
if [ "$EXITCODE" -ne 0 ]
then
  if [ "$EXITCODE" -ne 123 ]
  then
    echo -e "[$(date)]:\e[31m Exit code error from estimation.sh. Terminating\e[0m" >&2
    exit 1
  else
    echo -e "[$(date)]: \e[31mWARNING: Thread ran over time limit, will be re-processed in merging step...\e[0m" >&2
  fi
fi

# Merge estimation blocks
echo "[$(date)]: Call: mergeEstimation.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
mergeEstimation.sh -r "$refdir" -m "$matdir" -e "$tmpdir" -d "$estdir" -- "$mat_files" "$prefix" "$chr_num" "$total_part"

# Check if everything OK
if [ $? -ne 0 ] 
then
  echo -e "[$(date)]:\e[31m Nonzero exit code from mergeEstimation.sh Terminating...\e[0m" >&2
  exit -1
fi

# Run analysis on chunks
SECONDS=0
echo "[$(date)]: Call: singleMethAnalysis.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
seq "$total_part" | xargs -I {X} --max-proc "$threads" bash -c "timeout --signal=SIGINT '$time_limit'm singleMethAnalysis.sh -r '$refdir' -m '$matdir' -e '$estdir' -d '$tmpdir' --MC --ESI --MSI -- '$prefix' '$chr_num' '$total_part' {X}"

# Check if everything OK
EXITCODE="$?"
if [ "$EXITCODE" -ne 0 ]
then
  if [ "$EXITCODE" -ne 123 ]
  then
    echo -e "[$(date)]:\e[31m Nonzero exit code in singleMethAnalysis.sh. Terminating\e[0m" >&2
    exit 1
  else
    echo -e "[$(date)]: \e[31mWARNING: Thread ran over time limit, will be re-processed in merging step...\e[0m" >&2
  fi
fi

# Merge analysis blocks
echo "[$(date)]: Call: mergeSingleMethAnalysis.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
mergeSingleMethAnalysis.sh -r "$refdir" -e "$estdir" -a "$tmpdir" -d "$outdir" --MC --ESI --MSI -- "$prefix" "$chr_num" "$total_part"

# Check if everything OK
if [ $? -ne 0 ] 
then
  echo -e "[$(date)]:\e[31m Nonzero exit code in mergeSingleMethAnalysis.sh. Terminating...\e[0m" >&2
  exit -1
fi

# Time elapsed after succesful command
time_elapsed="$SECONDS"
echo "[$(date)]: Command succesful. Total time elapsed: $(( $time_elapsed / 3600)) h \
$(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
