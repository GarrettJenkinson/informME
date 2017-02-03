% informME: An information-theoretic pipeline for WGBS data
% Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                   objFnToMinimize.m                       %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 12/01/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function evaluates the objective function that must be minimized in 
% order to obtain maximum-likelihood estimates for the alpha, beta, and 
% gamma parameters of the Ising model. 
%
% USAGE:
%
% objFn = objFnToMinimize(theta)
%
% INPUT:
%
% theta
%       Parameter vector containing the alpha, beta, and gamma parameters 
%       of the Ising model.
%
% OUTPUT:
%
% objFn
%       The negative of the average "marginalized" log-likelihood of 
%       the input parameters given observations of the methylation 
%       state within a genomic region used for estimation. 
%

function objFn = objFnToMinimize(theta)

% Set global variables.

global CpGstart CpGend density Dist newMatrix 

% Compute model from parameters.

[An,Cn] = computeAnCn(density,Dist,theta);

% Compute average "marginalized" log-likelihood.

tempMat   = int32(newMatrix');
tempStart = int32(CpGstart);
tempEnd   = int32(CpGend);

aveLogLikelihood = computeAveLogLikelihood(An,Cn,tempMat,tempStart,tempEnd);

% Objective function to be minimized for maximum-likelihood estimation.

objFn = -aveLogLikelihood;