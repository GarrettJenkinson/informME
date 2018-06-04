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
matlab_function="mergeEstParams"

# Getopt command
TEMP="$(getopt -q -o hr:b:m:e:d:l: -l help,refdir:,bamdir:,matdir:,estdir:,outdir:,MATLICENSE: -n "$script_name.sh" -- "$@")"

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
outdir="$INTERDIR"
matdir="$INTERDIR"
estdir="$INTERDIR"
outdir="$INTERDIR"

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
    -m|--matdir)
      matdir="$2"
      shift 2
      ;;
    -e|--estdir)
      estdir="$2"
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
if [ "$#" -ne 4 ]; then
   echo -e "[$(date)]: \e[31mERROR: Command not valid. Check usage ...\e[0m" >&2
   cat "$script_absdir"/${script_name}_help.txt
   echo "[$(date)]: Terminating" >&2
   exit -1
fi
mat_files="$1"
prefix="$2"
chr_num="$3"
total_proc="$4"

# Check valid outdir
if [ -z "$outdir" ];then
   echo -e "[$(date)]: \e[31mERROR: Output directory is empty string ...\e[0m" >&2
   echo "[$(date)]: Terminating" >&2
   exit -1
fi

# Check valid directories
need_check=("$refdir" "$matdir" "$estdir")
informME_check_dirs "${need_check[@]}"

# Exit if previous function returned error
if [ $? -ne 0 ]
then
   echo "[$(date)]: Terminating" >&2
   exit -1
fi

# Output directory
mkdir -p "$outdir"
mkdir -p "${outdir}/chr${chr_num}"

# Set start time and message
SECONDS=0
echo "[$(date)]: Starting ..." 
echo "[$(date)]: Processing chromosome: ${chr_num}" 
echo "[$(date)]: Reference found in: ${refdir}" 
echo "[$(date)]: Matrices found in: ${matdir}" 
echo "[$(date)]: Estimations found in: ${estdir}" 

# Generate command and options
cmd="${matlab_function}('$mat_files','$prefix','$matdir','$refdir','$estdir','$chr_num','totalProcessors',$total_proc,'outdir','$outdir')"
options="-nodesktop -singleCompThread -nojvm -nosplash -nodisplay "

# Add license in case it is provided
if [ -n "${MATLICENSE}" ]
then
    options+="-c ${MATLICENSE}"
fi

# Start command
echo "[$(date)]: Starting command: ${cmd}"
matlab "${options}" -r "addpath('$matlab_library');try,${cmd},catch ME,fprintf(2,'Error identifier: %s. ',ME.identifier),fprintf(2,'Error message: %s\n',ME.message),exit(1),end,rmpath('$matlab_library'),exit(0)" 1>/dev/null

# Check if error in MATLAB, otherwise declare success
EXITCODE=$?
if [ $EXITCODE -ne 0 ]
then
   echo -e "[$(date)]: \e[31mERROR: MATLAB returned an error. Exiting ...\e[0m" >&2
   exit $EXITCODE
fi

# Time elapsed after succesful command
time_elapsed="$SECONDS"
echo "[$(date)]: Command succesful. Total time elapsed: $(( $time_elapsed / 3600)) h \
$(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
