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
%%%%%%%%%%%             Last Modified: 05/24/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the density of the n-th CpG site using the sorted 
% CpGlocation vector. 
%
% Usage:
% density = nDensity(CpGlocation,n,L)
%
% INPUT: 
%
% CpGlocation - Sorted vector of CpG locations on the current chromosome.
%
% n           - Index of CpG site.
%
% L           - Window size for density calculation.
%
% OUTPUT:
%
% density     - Computed density of n-th CpG site.
%

function density = nDensity(CpGlocation,n,L)

%
% compute upper and lower bounds of CpG sites which are included in 
% density calculation
%
UpperBound = CpGlocation(n)+floor(L/2);
LowerBound = CpGlocation(n)-floor(L/2);

%
% find the smallest and largest indices of the CpG sites between the
% upper and lower boundary (inclusive the boundary)
%
[lower_index,upper_index] = findSortedIndices(CpGlocation,LowerBound,...
                                              UpperBound);
                                   % CpGlocation must be sorted low to high

if isempty(lower_index)            % should not happen since current CpG 
                                   % always inside window 
    NumberOfCpGsBetweenBounds = 1; % current CpG always inside window
    disp('Warning: No CpGs found in density window.')
else
    NumberOfCpGsBetweenBounds = 1 + upper_index - lower_index; 
                        % +1 because inclusive at boundaries
end

density = double(NumberOfCpGsBetweenBounds)/double(L); 
                        % use double precision to represent density