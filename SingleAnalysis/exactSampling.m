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
%%%%%%%%                     exactSampling.m                       %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/04/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function draws an exact Monte Carlo sample from the Ising model 
% within a genomic region containing N CpG sites.
%
% USAGE:
%
% Xsamp = exactSampling(p1,transProbs)
%
% INPUTS:
%
% p1        
%            The probability Pr[X_1=1] of the first CpG site within the 
%            genomic region to be methylated.
%
% transProbs 
%            (N-1)x2 matrix whose elements are given by:  
%               transProbs(n,1) = Pr[X_{n+1}=0 | X_n=0]
%               transProbs(n,2) = Pr[X_{n+1}=0 | X_n=1]
%            where X_n is the methylation status of the n-th CpG 
%            site within the genomic region.
%
% OUTPUT:
%
% Xsamp
%           An Nx1 vector containg a sample drawn from the Ising model.
%

function Xsamp = exactSampling(p1,transProbs)

% Initialization.

N = size(transProbs,1)+1;
Xsamp = zeros(N,1);

% Draw an exact sample.

Xsamp(1) = rand <= p1;

for n = 1:(N-1)
    if Xsamp(n)==1
        Xsamp(n+1) = rand > transProbs(n,2);
    else
        Xsamp(n+1) = rand > transProbs(n,1);       
    end
end

