% informME: An information-theoretic pipeline for WGBS data
% Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu), 
% and Jordi Abante (jabante1@jhu.edu)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
% or see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                      compile.m                            %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 09/13/2017                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates MEX files from the C++ files located
% in the directory "MexSource". This requires configuration of
% a MATLAB compiler, as well as working instalations of the
% EIGEN (eigen.tuxfamily.org) and MPFR (http://www.mpfr.org)  
% packages. The generated binary MEX files are automatically
% placed in matlab_src.
%
% USAGE:
%
% compile(eigen_path,mpfr_path,mpfr_lib,cpp_library)
%
% INPUTS:
%
% eigen_path
%          Directory containing the files of the eigen_path instalation
%          (PATH/eigen_path).
%
% mpfr_path
%          Directory containing the "include" files of the MPFR
%          instalation (PATH/mpfr/include).
%
% mpfr_lib
%          Directory containing the "library" files of the MPFR
%          instalation (PATH/mpfr/lib).
%
% cpp_library
%          Directory containing the Matlab Mex C++ code files.

function compile(eigen_path,mpfr_path,mpfr_lib,cpp_library)

mex('-v','-largeArrayDims',['-I' mpfr_path],['-I' eigen_path],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp',strcat(cpp_library,'computeZ.cpp'))
mex('-v','-largeArrayDims',['-I' mpfr_path],['-I' eigen_path],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp',strcat(cpp_library,'computeZtilde.cpp'))
mex('-v','-largeArrayDims',['-I' mpfr_path],['-I' eigen_path],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp',strcat(cpp_library,'computeAveLogLikelihood.cpp'))
mex('-v','-largeArrayDims',['-I' mpfr_path],['-I' eigen_path],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp',strcat(cpp_library,'computeMCtransProbs.cpp'))
mex('-v','-largeArrayDims',['-I' mpfr_path],['-I' eigen_path],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp',strcat(cpp_library,'calcMargProb.cpp'))

