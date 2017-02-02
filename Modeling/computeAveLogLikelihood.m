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
% This function computes the average log likelihood of a set of
% observations for a given Ising model.
%
% avelogLikelihood = computeAveLogLikelihood(An,Cn,dataMatrix, CpGstart, CpGend)
%
% Inputs:
%
% An is a Nx1 vector of a_n parameters for n=1,...,N
%
% Cn is a (N-1)x1 vector of C_n parameters for n=1,...,N-1
%
% dataMatrix
%           A matrix with each row corresponding to a read and each column
%           corresponding to a CpG site. It takes values of -1, 0, 1
%           which indicate that a given read does not observe a CpG site,
%           observes a lack of methylation at a CpG site, or observes
%           methylation at a CpG site, respectively.
%
% CpGstart
%           A vector with as many elements as rows of newMatrix. Each
%           element contains the index of the first observation in the
%           corresponding row of newMatrix.
%
% CpGend
%           A vector with as many elements as rows of newMatrix. Each
%           element contains the index of the last observation in the
%           corresponding row of newMatrix.
%
% Outputs:
function avelogLikelihood = computeAveLogLikelihood(An,Cn,dataMatrix, CpGstart, CpGend)

dataMatrix = dataMatrix'; % get matrix with columns each read, convert back to rows

d = 100; % number of digits for variable precision arithmetic

%
% Compute Partition Function
%

[logZ1,logZ0,logZ]      = computeZ(An,Cn);
[logZ1tilde,logZ0tilde] = computeZtilde(An,Cn); 

%
%  Compute average log-likelihood
%

avelogLikelihood = sym(0);

K = length(CpGstart);
for k=1:K
    y_d              = dataMatrix(k,CpGstart(k):CpGend(k));
    
    logMargProb      = calcMargProb(CpGstart(k),CpGend(k)-CpGstart(k),y_d,...
                                    logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,...
                                    An,Cn);
    
    avelogLikelihood = vpa(avelogLikelihood + logMargProb,d);
end

avelogLikelihood = double(vpa(avelogLikelihood/K,d));