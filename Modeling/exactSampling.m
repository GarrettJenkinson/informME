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
% A function that draws an exact Monte Carlo sample from the Ising model.
%
% Xsamp = exactSampling(p1,transProbs)
%
% INPUTS:
% p1        
%           Scalar whose value is Pr[X_1=1]
%
% transProbs 
%           (N-1)-by-2 Matrix an such that: 
%           transProbs(n,1) is Pr[x_{n+1}=0|x_{n}=0]
%           transProbs(n,2) is Pr[x_{n+1}=0|x_{n}=1]
%
% OUTPUTS:
%
% Xsamp
%           A Nx1 vector that is a Monte Carlo sampling from Pr[X]

function Xsamp = exactSampling(p1,transProbs)


%
% Initialize
%

N = size(transProbs,1)+1;

Xsamp = zeros(N,1);

%
% Compute an exact sample from the inhomogenous markov chain
%
Xsamp(1) = rand<=p1; % p1 is prob Xsamp(N)=1, 
                     % rand<=p1 occurs with probability p1

for n=1:(N-1)
    if Xsamp(n)==1 %X_{n}=1
        
        Xsamp(n+1) = rand>transProbs(n,2); %transProbs(n,2) is p(x_{n+1}=0|x_{n}=1)
                                           % Xsamp(n+1) is 1 with probability 1-transProbs(n,2)=p(x_{n+1}=1|x_{n}=1)
                                            
        
    else          %X_{n}=0
        
        Xsamp(n+1) = rand>transProbs(n,1);%transProbs(n,1) is p(x_{n+1}=0|x_{n}=0)
                                          % Xsamp(n+1) is 1 with probability 1-transProbs(n,2)=p(x_{n+1}=1|x_{n}=0)
        
    end
end