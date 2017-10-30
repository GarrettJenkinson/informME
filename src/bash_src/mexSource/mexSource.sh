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

# Find matlab directory 3 levels above
aux="$(dirname "$(dirname "$script_absdir")")"
matlab_library="${aux}/matlab_src/"
cpp_library="${aux}/cpp_src/"
matlab_function="compile"

# Print help if no arguments
if [ $# -eq 0 ]
    then
        cat "$script_absdir/${script_name}_help.txt"
        exit 1
fi

# Getopt command
TEMP="$(getopt -o hl: -l help,MATLICENSE: -n "$script_name.sh" -- "$@")"

if [ $? -ne 0 ] 
then
  echo "Terminating..." >&2
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
    -l|--MATLICENSE)
      MATLICENSE="$2"
      shift 2
      ;;  
    --) 
      shift
      break
      ;;  
    *)  
      echo "$script_name.sh:Internal error!"
      exit -1
      ;;  
  esac
done

# Get inputs
eigen_path="$1"
mpfr_inc_path="$2"
mpfr_lib_path="$3"

# Set start time and message
SECONDS=0
echo "[$(date)]: Checking ..." 

# Generate command and options
cmd="${matlab_function}('$eigen_path','$mpfr_inc_path','$mpfr_lib_path','$cpp_library')"
options="-nodesktop -singleCompThread -nojvm -nosplash -nodisplay "

# Add license in case it is provided
if [ -n "${MATLICENSE}" ]
then
    options+="-c ${MATLICENSE}"
fi

# Start command
echo "[$(date)]: Starting command: ${cmd}" 
matlab "${options}" -r "addpath('$matlab_library');${cmd};rmpath('$matlab_library');exit"

# Move .mex files to matlab_src
mv *.mex* "$matlab_library"

# Check if error in MATLAB, otherwise declare success
if [ $? -ne 0 ]
then
   echo "Error in MATLAB"
   echo "$cmd"
   exit $EXITCODE
fi

# Time elapsed after succesful command
time_elapsed="$SECONDS"
echo "[$(date)]: Command succesful. Total time elapsed: $(( $time_elapsed / 3600)) h \
$(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
