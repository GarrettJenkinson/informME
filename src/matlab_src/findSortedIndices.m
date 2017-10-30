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
%%%%%%%%                  findSortedIndices.m                      %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/30/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function uses binary search to rapidly find the smallest and the 
% largest indices of the elements of a vector of sorted values (low to 
% high) that are between a lower and an upper bound.
%
% USAGE:
%
% [lower_index,upper_index] = findSortedIndices(xvec,LowerBound,UpperBound)
%
% INPUTS:
%
% xvec
%            Vector of sorted values (low to high).
%
% LowerBound
%            Lower bound on xvec values.
%
% UpperBound
%            Upper bound on xvec values.
%
% OUTPUTS:
%
% lower_index
%            Smallest index such that LowerBound <= xvec(index) <= UpperBound.
%
% upper_index
%            Largest index such that LowerBound <= xvec(index) <= UpperBound.
%

function [lower_index,upper_index] = findSortedIndices(xvec,LowerBound,UpperBound)

if LowerBound > xvec(end) || UpperBound < xvec(1) || UpperBound < LowerBound
    % No indices satify bounding condition. 
    lower_index = [];
    upper_index = [];
    return;  
end

lower_index_a = int64(1);         
lower_index_b = int64(length(xvec)); 
                         % xvec(lower_index_b) satisfies lowerbound.
upper_index_a = int64(1);         
                         % xvec(upper_index_a) satisfies upperbound.
upper_index_b = int64(length(xvec)); 

% The following loop increases lower_index_a and upper_index_a and 
% decreases lower_index_b and upper_index_b until they differ by at 
% most 1. Because one of these index variables always satisfies the 
% appropriate bound, this means that the loop will terminate with 
% either lower_index_a or lower_index_b having the minimum possible 
% index that satifies the lower bound, or upper_index_a or upper_index_b 
% having the largest possible index that satisfies the upper bound. 

while (lower_index_a+1 < lower_index_b) || (upper_index_a+1 < upper_index_b)
    lw = (lower_index_a+lower_index_b)/int64(2); % Split the lower index.
    
    if xvec(lw) >= LowerBound
        lower_index_b=lw; % Decrease lower_index_b (whose xvec value 
                          % remains >= to lower bound).   
    else
        lower_index_a=lw; % Increase lower_index_a.
        if (lw>upper_index_a) && (lw<upper_index_b)
            upper_index_a=lw; % Increase upper_index_a (whose xve  value 
                              % remains < to lower bound and thus to 
                              % the upper bound).
        end
    end

    up = (upper_index_a+upper_index_b)/int64(2); % split the lower index
    if xvec(up) <= UpperBound
        upper_index_a=up;     % Increase upper_index_a (whose xvec value  
                              % remains <= to upper bound). 
    else
        upper_index_b=up;     % Decrease upper_index_b.
        if (up<lower_index_b) && (up>lower_index_a)
            lower_index_b=up; % Decrease lower_index_b (whose xvec value 
                              % remains > to upper bound and thus to 
                              % the lower bound).
        end
    end
end

% Choose the final answer.

if xvec(lower_index_a) >= LowerBound
    lower_index = lower_index_a;
else
    lower_index = lower_index_b;
end
if xvec(upper_index_b) <= UpperBound
    upper_index = upper_index_b;
else
    upper_index = upper_index_a;
end

if lower_index > upper_index
    lower_index=[];
    upper_index=[];
end
