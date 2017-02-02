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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Statistical Model for DNA Methylation Patterns    %%%%%%%%%%%%
%%%%%%%%%%%  Code by: Garrett Jenkinson                        %%%%%%%%%%%%
%%%%%%%%%%%             Last Modified: 05/15/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function computes the objective function that must be minimized in
% optimization in order to produce the maximum likelihood statistical
% estimate of the parameters. 
%
% objFn = objFnToMinimize(theta);
%
% It takes as input:
%
% theta
%       The parameter vector that produces the Ising model for the region.
%
% And gives as output:
%
% objFn
%       The negative average log likelihood of the input theta under the
%       observed methylation sample. Choosing a theta that minimizes this
%       value will produce a maximum likelihood parameter estimate.

function objFn = objFnToMinimize(theta)

% set global vars
global CpGstart CpGend density Dist newMatrix 

%
% Compute model from parameters
%
[invOmega_An,invOmega_Cnm] = computeAnCnm(density,Dist,theta);

%
% Compute average log likelihood (function to maximize)
%
tempMat = int32(newMatrix');
tempStart = int32(CpGstart);
tempEnd = int32(CpGend);

aveLogLikelihood = computeAveLogLikelihood(invOmega_An,invOmega_Cnm,tempMat, tempStart, tempEnd);

%
% Objective function to minimize is the negative of the function to
% maximize
%
objFn = -aveLogLikelihood;
