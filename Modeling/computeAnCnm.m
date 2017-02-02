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
% This function takes in the density and distances between CpG sites in a 
% region, along with the theta vector and computes the a_n and 
% c_nm parameters for the Ising model of the region.
%
% [An,Cnm] = computeAnCnm(density,Dist,theta)
%
% This takes as inputs:
%
% theta   
%           ( \alpha, \beta, \gamma, a_1/\Omega,a_N/\Omega)
%           Contains the 5 parameters of the model           
% density  
%           \rho_n for n\in\sN
%           Contains the local CpG density for each CpG site
% Dist
%           d(n,n') for n \in {1,2,...,N-1}, and n'=n+1
%           Contains the basepair distances between each adjacent CpG site
%
% And provides outputs:
%
% An
%           The vector of length N of a_n parameters
%
% Cnm
%           The vector of length N-1 of c_nm parameters

function [An,Cnm] = computeAnCnm(density,Dist,theta)


%
% a_n = \alpha+\beta*density; n=2,...,N-1
% a_n = a_n; n=1,N
%
An      = theta(1)+(theta(2)*density);
An(1)   = theta(4);
An(end) = theta(5);

%
% Cnm = gamma./(Dist);
%
Cnm = theta(3)./double(Dist);