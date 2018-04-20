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

#
# Specify what versions of dependencies to download/install
#

# tar.bz2 sources
GMPSRC="http://mirror.team-cymru.org/gnu/gmp/gmp-6.1.2.tar.bz2"

# tar.gz sources
MPFRSRC="http://www.mpfr.org/mpfr-3.1.5/mpfr-3.1.5.tar.gz"
EIGSRC="http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz"
MPREALSRC="https://bitbucket.org/advanpix/mpreal/get/mpfrc++-3.6.5.tar.gz"


#
# Make sure installation script dependencies are met
#

clear 

printf "InformME  Copyright (C) 2017  Garrett Jenkinson and Jordi Abante\n    This program comes with ABSOLUTELY NO WARRANTY; for details type \'cat GPL3_LICENSE.txt\'.\n    This is free software, and you are welcome to redistribute it\n    under certain conditions; type \'cat GPL3_LICENSE.txt\' for details.\n\n"

printf "Hello $USER, the informME installation is starting...\n"

STARTDIR=$(pwd)
INFORMMEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

sleep 2

printf "Before running this script, you must have setup matlab's mex compiler.\n\nSee here for details: \n https://www.mathworks.com/help/matlab/ref/mex.html \n\n"

read -p "Is Matlab mex configured for this machine (y/n)? " -n 1 -r

if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    printf "\n\nExiting. Rerun script after configuring mex.\n"
    exit 1
fi

clear

printf "\n\nWe will begin by installing any dependencies that are not yet installed. If you are on a cluster, this script needs to be run on a compute node and NOT the login node.\n\n"

read -p "Are you running on a node that can be used for computation (y/n)? " -n 1 -r

if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    printf "\n\nExiting. Rerun script from compute node.\n"
    exit 2
fi

clear


#
# Get info or install GMP and MPRF
#
printf "Throughout this process, you may be asked to provide directories. Note that tab-completion is enabled in this process to make the process easier. \n\n"

read -p "Are GMP and MPFR installed already (y/n)? " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then
    printf "\nGreat! Please enter the path to the GMP install (and hit enter): \n"
    read -e GMPDIR

    GMPDIR=$(readlink -f "`eval echo ${GMPDIR//>}`")

    if [ ! -r "$GMPDIR" ] 
    then
        printf "$GMPDIR is not readable\n"
        exit 3
    fi
    
    printf "\nPlease enter the path to the MPFR install (and hit enter): \n"
    read -e MPFRDIR

    MPFRDIR=$(readlink -f "`eval echo ${MPFRDIR//>}`")

    if [ ! -r "$MPFRDIR" ] 
    then
        printf "$MPFRDIR is not readable\n"
        exit 4
    fi
    clear
else
    printf "\n No problem, we will try to install them.\n\nPlease enter folder should we install it in (and hit enter):\n"
    read -e GMPDIR
    
    GMPDIR=$(readlink -f "`eval echo ${GMPDIR//>}`")

    if [ ! -w "$GMPDIR" ] 
    then
        printf "$GMPDIR is not writable\n"
        exit 5
    fi
    
    cd $GMPDIR 
    mkdir informMEdeps
    GMPDIR=$GMPDIR/informMEdeps
    MPFRDIR=$GMPDIR
    cd $GMPDIR
    mkdir bin
    mkdir include
    mkdir lib
    mkdir share
    mkdir src 
    cd src
    wget $GMPSRC
    wget $MPFRSRC
    mkdir $GMPDIR/src/gmp
    tar -xjf *gmp*tar* -C $GMPDIR/src/gmp --strip-components=1
    mkdir $GMPDIR/src/mpfr
    tar -xzf *mpfr*tar* -C $GMPDIR/src/mpfr --strip-components=1
    rm *tar*
    cd gmp
    ./configure --prefix=$GMPDIR
    make
    make install
    cd ../mpfr
    ./configure --prefix=$MPFRDIR --with-gmp=$GMPDIR
    make
    make install
    
    clear
fi


#
# Get info or install Eigen and mpreal.h
#


read -p "Are Eigen and mpreal.h installed already (y/n)? " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then
    printf "\nGreat! Please enter the path to the eigen install (and hit enter): \n"
    read -e EIGDIR

    EIGDIR=$(readlink -f "`eval echo ${EIGDIR//>}`")

    if [ ! -r "$EIGDIR" ] 
    then
        printf "$EIGDIR is not readable\n"
        exit 6
    fi
    
    clear
else
    printf "\n No problem, we will try to install them.\n\nPlease enter the folder should we install it in (and hit enter):\n"
    read -e EIGDIR
    
    EIGDIR=$(readlink -f "`eval echo ${EIGDIR//>}`")    

    if [ ! -w "$EIGDIR" ]
    then
        printf "$EIGDIR is not writable\n"
        exit 7
    fi

    cd $EIGDIR
    if [ ! -w "$EIGDIR/informMEdeps" ]
    then
        mkdir informMEdeps
    fi
    cd informMEdeps
    wget --no-check-certificate $EIGSRC    
    EIGDIR=$EIGDIR/informMEdeps/eigen
    mkdir $EIGDIR
    tar -xzf *tar.gz -C $EIGDIR --strip-components=1
    cd $EIGDIR
    cp unsupported/Eigen/MPRealSupport Eigen/MPRealSupport
    wget --no-check-certificate $MPREALSRC 
    mkdir $EIGDIR/mpreal
    tar -xzf *mp*tar.gz -C $EIGDIR/mpreal --strip-components=1
    cp $EIGDIR/mpreal/mpreal.h $EIGDIR/mpreal.h
    clear
fi


cd $STARTDIR

#
# Make changes to user's .bashrc
#

printf '\ninformME would like to make the following additions to your .bashrc file\n\n\n#Added by informME installation\nif [ -z "${LD_LIBRARY_PATH}" ]; then\n\texport LD_LIBRARY_PATH="'
printf "$MPFRDIR/lib"
printf '"\nelse\n\texport LD_LIBRARY_PATH="'
printf "$MPFRDIR/lib"
printf ':${LD_LIBRARY_PATH}"\nfi\nexport PATH="$PATH:'
printf "${INFORMMEDIR}"
printf '/bin"\n#End added by informME\n\n'

read -p "Do you want these lines appended to your ~/.bashrc file (y/n)? " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then
    printf '\n\n\n#Added by informME installation\nif [ -z "${LD_LIBRARY_PATH}" ]; then\n\texport LD_LIBRARY_PATH="' >> ~/.bashrc
    printf "$MPFRDIR/lib" >> ~/.bashrc
    printf '"\nelse\n\texport LD_LIBRARY_PATH="' >> ~/.bashrc
    printf "$MPFRDIR/lib" >> ~/.bashrc
    printf ':${LD_LIBRARY_PATH}"\nfi\nexport PATH="$PATH:' >> ~/.bashrc
    printf "${INFORMMEDIR}" >> ~/.bashrc
    printf '/bin"\n#End added by informME\n\n' >> ~/.bashrc
fi

clear


#
# Create informME hidden folder and config file 
#

printf "informME can create a configuration file ~/.informME/informME.config that will save default input/output folders and dependency information so that you don't need to specify this information every time you use the tool.\n\n"

printf "Later, if you specify a folder location as an argument to an informME script, the script will ignore the config file settings. Also the config file can be edited/updated at any time.\n\n"

read -p "Do you want to create this config file (y/n)? " -n 1 -r

if [[ $REPLY =~ ^[Yy]$ ]]
then   
    printf "\n\n Editing config file...\n\n We will need some input from you now. If you do not know the location or do not have a default location in mind, then just specify ~ for the directory, and later you will have to edit the config file or pass the correct information to the appropriate informME script.\n\n"
    mkdir ~/.informME
    cd ~/.informME
    
    # note that if a config file exists, we do not erase but
    # simply append to that existing, which overwrites the old
    # settings, but allows someone to see the old settings
    # and revert as desired.
    printf "# install.sh generated config file\n" >> informME.config
    printf "# generated on: $(date +%m-%d-%Y)\n\n" >> informME.config
    printf "MPFRDIR=${MPFRDIR}\n" >> informME.config
    printf "EIGDIR=${EIGDIR}\n" >> informME.config
    
    # get scratch disk information
    printf "\n\nA scratch disk is a location where we can save temporary files to disk and perform frequent read/write operations quickly. The data on the scratch disk will not be stored long term, but it does need to persist beyond each individual job's lifetime.\n\n"    
    
    printf "\nPlease enter the scratch disk folder path (and hit enter): \n"
    read -e SCRATCHDIR

    SCRATCHDIR=$(readlink -f "`eval echo ${SCRATCHDIR//>}`")

    if [ ! -w "$SCRATCHDIR" ]
    then
        printf "$SCRATCHDIR is not writeable\n"
        exit 8
    fi
    printf "SCRATCHDIR=${SCRATCHDIR}\n" >> informME.config
    
    # get output fasta analysis location
    printf "\n\nThe output folder where informME will store results of analyzing a reference genome. This folder will have subdirectories for each reference genome. So for example, if the folder you specify is \n/path/to/refGenomes\n then when you use informME to analyze hg19.fasta informME store its results from analyzing the reference genome in\n/path/to/refGenomes/hg19/\nwhich will able to be re-used for all samples on that reference genome (i.e., each reference genome is only analyzed one time). \n\n"

    printf "\nPlease enter the folder path to store these results (and hit enter): \n"
    read -e REFGENEDIR

    REFGENEDIR=$(readlink -f "`eval echo ${REFGENEDIR//>}`")

    if [ ! -w "$REFGENEDIR" ]
    then
        printf "$REFGENEDIR is not writeable\n"
        exit 9
    fi
    printf "REFGENEDIR=${REFGENEDIR}\n" >> informME.config
    
    # get input indexed BAM folder 
    printf "\n\nThe input folder where informME will find the bam files and corresponding index bai files from the WGBS experiments.\n\n"

    printf "\nPlease enter the folder path that will provide these bam files (and hit enter): \n"
    read -e BAMDIR

    BAMDIR=$(readlink -f "`eval echo ${BAMDIR//>}`")

    if [ ! -r "$BAMDIR" ]
    then
        printf "$BAMDIR is not readable\n"
        exit 10
    fi
    printf "BAMDIR=${BAMDIR}\n" >> informME.config
    

    # get output intermediate .mat folder 
    printf "\n\nThe output folder where informME will store its intermediate .mat file results. These files are useful in that they store the matrices processed from bam files or the statistical models built from these data. By storing them, it is not necessary to recompute things when comparing against a new sample or building a new pooled model, etc. These files can take up a lot of space and need to persist for as long as you might reanalyze a given dataset.\n\n"

    printf "\nPlease enter the folder path that will store these files (and hit enter): \n"
    read -e INTERDIR

    INTERDIR=$(readlink -f "`eval echo ${INTERDIR//>}`")

    if [ ! -w "$INTERDIR" ]
    then
        printf "$INTERDIR is not writeable\n"
        exit 11
    fi
    printf "INTERDIR=${INTERDIR}\n" >> informME.config


    # get final output folder 
    printf "\n\nThe output folder where informME will store its final results.\n\n"

    printf "\nPlease enter the folder path that will store these files (and hit enter): \n"
    read -e FINALDIR

    FINALDIR=$(readlink -f "`eval echo ${FINALDIR//>}`")

    if [ ! -w "$FINALDIR" ]
    then
        printf "$FINALDIR is not writeable\n"
        exit 12
    fi
    printf "FINALDIR=${FINALDIR}\n" >> informME.config

    
    # get final output folder 
    printf "\n\nThe matlab license file (.lic), which allows Matlab to run and use the Bioinformatics and Symbolic Math toolboxes. It is better to specify this file when it is known. \n\n"
    
    read -p "Do you want to specify the matlab license file (y/n)? " -n 1 -r

    if [[ $REPLY =~ ^[Yy]$ ]]
    then

        printf "\nPlease enter the full path and file name (e.g., /path/to/file.lic) of the matlab license (and hit enter): \n"
        read -e MATLICENSE

        MATLICENSE=$(readlink -f "`eval echo ${MATLICENSE//>}`")

        if [ ! -r "$MATLICENSE" ]
        then
            printf "$MATLICENSE is not readable\n"
            exit 13
        fi
        printf "MATLICENSE=${MATLICENSE}\n" >> informME.config

    fi

    clear
fi


#
# Compile the mex code
#


if [ -r "$MATLICENSE" ]
then
    ${INFORMMEDIR}/bin/mexSource.sh -l ${MATLICENSE}  -- ${EIGDIR} ${MPFRDIR}/include ${MPFRDIR}/lib
else
    ${INFORMMEDIR}/bin/mexSource.sh -- ${EIGDIR} ${MPFRDIR}/include ${MPFRDIR}/lib
fi
cd $STARTDIR 


printf "\nInformME is now installed. Please see the readme for more information on using the tool.\n\n"

