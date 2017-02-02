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
libpath = '/home/student/wjenkin6/MexCode/mpfr/lib';
headpath2 = '/home/student/wjenkin6/MexCode/mpfr/include';
headpath1 = '/home/student/wjenkin6/MexCode/eigen';


mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeZ.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeZtilde.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeAveLogLikelihood.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','computeMCtransProbs.cpp')

mex('-v','-largeArrayDims',['-I' headpath1],['-I' headpath2],['-L' libpath],'-lmpfr','-lgmp','calcMargProb.cpp')
