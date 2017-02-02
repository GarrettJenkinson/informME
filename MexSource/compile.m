%   informME: An information-theoretic pipeline for WGBS data
%   Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
%   or see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                      compile.m                            %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/29/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates MEX files from the C++ files located
% in the directory "MexSource". This requires configuration of
% a MATLAB compiler, as well as working instalations of the
% EIGEN (eigen.tuxfamily.org) and MPFR (http://www.mpfr.org)  
% packages. The generated binary MEX files are automatically
% placed in the "Modeling" and "SingleAnalysis" directories of 
% informME.
%
% USAGE:
%
% compile(eigen,mpfr_incl,mpfr_lib,eigen)
%
% INPUTS:
%
% eigen
%          Directory containing the files of the EIGEN instalation
%          (PATH/eigen).
%
% mpfr_inc
%          Directory containing the "include" files of the MPFR
%          instalation (PATH/mpfr/include).
%
% mpfr_lib
%          Directory containing the "library" files of the MPFR
%          instalation (PATH/mpfr/lib).
%

function compile(eigen,mpfr_incl,mpfr_lib)

mex('-v','-largeArrayDims',['-I' mpfr_incl],['-I' eigen],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp','computeZ.cpp')
mex('-v','-largeArrayDims',['-I' mpfr_incl],['-I' eigen],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp','computeZtilde.cpp')
mex('-v','-largeArrayDims',['-I' mpfr_incl],['-I' eigen],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp','computeAveLogLikelihood.cpp')
mex('-v','-largeArrayDims',['-I' mpfr_incl],['-I' eigen],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp','computeMCtransProbs.cpp')
mex('-v','-largeArrayDims',['-I' mpfr_incl],['-I' eigen],...
    ['-L' mpfr_lib],'-lmpfr','-lgmp','calcMargProb.cpp')

% Move the generated MEX files to correct directories.

status1 = system('cp *.mex* ../Modeling/');

if status1 ~= 0
    disp('Error moving MEX binary files to Modeling directory');
end

status2 = system('mv *.mex* ../SingleAnalysis/');

if status2 ~= 0
    disp('Error moving MEX binary files to SingleAnalysis directory');
end

end