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
%%%%%%%%                    computeAnCn.m                          %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 12/01/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the a_n and c_n parameters of the Ising model 
% within a genomic region by using the CpG densities, the CpG distances, 
% and the estimated alpha, beta and gamma parameters of the model. 
%
% USAGE:
%
% [An,Cn] = computeAnCn(density,distance,theta)
%
% INPUTS:
%
% theta   
%           Vector containing the estimated values of the five alpha, 
%           beta, and gamma parameters of the Ising model.
%
% density  
%           Vector containing the CpG densities of the CpG sites within 
%           the genomic region. 
%          
% distance
%           Vector containing the CpG distances associated with the 
%           CpG sites within the genomic region. 
%
% OUTPUTS:
%
% An
%           Vector containing the values of the a_n parameters.
%
% Cn
%           Vector containing the values of the c_n parameters.
%

function [An,Cn] = computeAnCn(density,distance,theta)

An      = theta(1)+(theta(2)*density);
An(1)   = theta(4);
An(end) = theta(5);

Cn = theta(3)./double(distance);
