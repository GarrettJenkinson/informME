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
   source informME_source_config
fi

# Find matlab directory 3 levels above
aux="$(dirname "$(dirname "$(dirname "$script_absdir")")")"
matlab_library="${aux}/matlab_src/"
matlab_function="$script_name"

# Getopt command
TEMP="$(getopt -q -o hr:b:d:t:c:p:l:q: -l help,refdir:,bamdir:,tmpdir:,outdir:,trim:,chr_string:,paired_ends:,MATLICENSE:,threads:,time_limit:,total_part: -n "$script_name.sh" -- "$@")"

if [ $? -ne 0 ]
then
  echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
  cat "$script_absdir"/${script_name}_help.txt
  echo "[$(date)]: Terminating" >&2
  exit -1
fi

eval set -- "$TEMP"

# Defaults
refdir="$REFGENEDIR"
bamdir="$BAMDIR"
tmpdir="$SCRATCHDIR"
outdir="$INTERDIR"
trim=0
chr_string=1
paired_ends=1
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
    -b|--bamdir)
      bamdir="$2"
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
    -t|--trim)
      trim="$2"
      shift 2
      ;;  
    -c|--chr_string)
      chr_string="$2"
      if ([ "$chr_string" -ne "0" ] && [ "$chr_string" -ne "1" ])
      then 
        echo -e "[$(date)]: \e[31mNot a valid choice of -c option, must be either 0 or 1. Terminating...\e[0m" >&2
        exit -1
      fi
      shift 2
      ;;  
    -p|--paired_ends)
      paired_ends="$2"
      if ([ "$paired_ends" -ne "0" ] && [ "$paired_ends" -ne "1" ])
      then
        echo -e "[$(date)]: \e[31mNot a valid choice of -p option, must be either 0 or 1. Terminating...\e[0m" >&2
        exit -1
      fi
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
bam_file="$1"
chr_num="$2"

# Get BAM file
bam_prefix="$(basename "$bam_file" .bam)"

# Check valid outdir
if [ -z "$outdir" ];then
   echo -e "[$(date)]: \e[31mERROR: Output directory is empty string ...\e[0m" >&2
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
mkdir -p "$outdir"
mkdir -p "${outdir}/chr${chr_num}"

# Check that input files exist
bamdir="$(readlink -f "`eval echo ${bamdir//>}`")/"
if [ ! -r "${bamdir}${bam_prefix}.bam" ]
then
  echo -e "[$(date)]: \e[31mInput bam file:\e[0m" >&2
  echo -e "[$(date)]:\e[31m ${bamdir}${bam_prefix}.bam\e[0m" >&2
  echo -e "[$(date)]:\e[31m is not readable. Terminating...\e[0m" >&2
  exit 1
fi
if [ ! -r "${bamdir}${bam_prefix}.bam.bai" ]
then
  echo -e "[$(date)]:\e[31m Input bam index file:\e[0m" >&2
  echo -e "[$(date)]:\e[31m ${bamdir}${bam_prefix}.bam.bai\e[0m" >&2
  echo -e "[$(date)]:\e[31m is not readable. Terminating...\e[0m" >&2
  exit 1
fi

# Check that samtools is available
if [ ! -x "`which samtools`" ]
then
  echo -e "[$(date)]:\e[31m Samtools not installed:\e[0m" >&2
  echo -e "[$(date)]:\e[31m samtools is not executable on system PATH\e[0m" >&2
  echo -e "[$(date)]:\e[31m Terminating...\e[0m" >&2
  exit 1
fi

#check if final output already exists
if [ -r "${outdir}/chr${chr_num}/${bam_prefix}_matrices.mat" ]
then
  echo -e "[$(date)]:\e[31m Warning: Final output file:\e[0m" >&2
  echo -e "[$(date)]:\e[31m ${outdir}/chr${chr_num}/${bam_prefix}_matrices.mat\e[0m" >&2
  echo -e "[$(date)]:\e[31m already exists. Delete file and rerun if you wish to overwrite.\e[0m" >&2
  echo -e "[$(date)]:\e[31m Terminating...\e[0m" >&2
  exit 0
fi

# Get matrices via matrixFromBam.sh
SECONDS=0
echo "[$(date)]: Call: matrixFromBam.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
seq "$total_part" | xargs -I {X} --max-proc "$threads" bash -c "timeout --signal=SIGINT '$time_limit'm matrixFromBam.sh -r '$refdir' -b '$bamdir' -d '$tmpdir' -c '$chr_string' -p '$paired_ends' -t '$trim' -- '$bam_prefix' '$chr_num' '$total_part' {X}"
 
# Check if everything OK or the job was interrupted due to excessive length
EXITCODE="$?"
if [ "$EXITCODE" -ne 0 ]
then
  if [ "$EXITCODE" -ne 123 ]
  then 
    echo -e "[$(date)]: \e[31m Error in exit code matrixFromBam.sh. Terminating\e[0m" >&2
    exit 1
  else
    echo -e "[$(date)]: \e[31mWARNING: Thread ran over time limit, will be re-processed in merging step...\e[0m" >&2
  fi
fi

# Merge matrices via mergeMatrices.sh
echo "[$(date)]: Call: mergeMatrices.sh ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
mergeMatrices.sh -r "$refdir" -b "$bamdir" -m "$tmpdir" -d "$outdir" -c "$chr_string" -p "$paired_ends" -t "$trim" -- "$bam_file" "$chr_num" "$total_part"

# Check if everything OK
if [ $? -ne 0 ] 
then
  echo -e "[$(date)]:\e[31m Error in exit code at mergeMatrices.sh. Terminating\e[0m" >&2
  exit -1
fi

# Time elapsed after succesful command
time_elapsed="$SECONDS"
echo "[$(date)]: Command succesful. Total time elapsed: $(( $time_elapsed / 3600)) h \
$(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
