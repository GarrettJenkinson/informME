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
% This function computes the E[S(X)] as well as the nearest neighbor 
% correlation coefficients for the Ising model.
%
% [ExpectedValSofX,NNcorr] = expectValSXexact(margProbs,transProbs,density,Dist)
%
% INPUTS:
%
% margProbs         A Nx1 vector of the marginal probabilities of methylation.
%
% transProbs        (N-1)-by-2 Matrix an such that: 
%                   transProbs(n,1) is Pr[x_{n+1}=0|x_{n}=0]
%                   transProbs(n,2) is Pr[x_{n+1}=0|x_{n}=1]
%
% Dist
%                   A (N-1)x1 vector of distances to the next CpG site for 
%                   each CpG site
%
% density           A Nx1 vector of densities for each CpG site
%
% OUTPUTS:
%
% ExpectedValSofX   A 5x1 vector for E[S(X)].
%
% NNcorr            A (N-1)x1 vector of nearest neighbor correlation coeffs


function [ExpectedValSofX,NNcorr] = expectValSXexact(margProbs,transProbs,density,Dist)

N=length(density);
ExpectedValSofX = zeros(5,1);
Dist = double(Dist);

%
% Compute Pearson Correlation Coeff
%
stdDevs = sqrt(margProbs.*(1-margProbs));

% E[(X_n-\mu_n)(X_{n+1}-\mu_{n+1})]/\sigma_n\sigma_{n+1}
NNcorr = ( margProbs(1:(N-1)).*margProbs(2:N).*(1-margProbs(1:(N-1))).*transProbs(1:(N-1),1)... % X_n=0 X_{n+1}=0
          - margProbs(1:(N-1)).*(1-margProbs(2:N)).*(1-margProbs(1:(N-1))).*(1-transProbs(1:(N-1),1))... % X_n=0 X_{n+1}=1
          - (1-margProbs(1:(N-1))).*margProbs(2:N).*margProbs(1:(N-1)).*transProbs(1:(N-1),2)  ... % X_n=1 X_{n+1}=0
          + (1-margProbs(1:(N-1))).*(1-margProbs(2:N)).*margProbs(1:(N-1)).*(1-transProbs(1:(N-1),2))  )... % X_n=1 X_{n+1}=1
        ./(stdDevs(1:(N-1)).*stdDevs(2:N));  
%
% Compute Expected values
%    

ExpectedValSofX(1)=sum(2*margProbs(2:(end-1))-1);%sum(2*Xsamp(2:(end-1))-1);
ExpectedValSofX(2)=sum( density(2:(end-1)).*(2*margProbs(2:(end-1))-1) );%density(2:(end-1))'*(2*Xsamp(2:(end-1))-1);
ExpectedValSofX(3)=sum( (1-margProbs(1:(N-1))).*transProbs(1:(N-1),1)./Dist... % X_n=0 X_{n+1}=0
                      - (1-margProbs(1:(N-1))).*(1-transProbs(1:(N-1),1))./Dist... % X_n=0 X_{n+1}=1
                      - margProbs(1:(N-1)).*transProbs(1:(N-1),2)./Dist  ... % X_n=1 X_{n+1}=0
                      + margProbs(1:(N-1)).*(1-transProbs(1:(N-1),2))./Dist  ); % X_n=1 X_{n+1}=1                        
ExpectedValSofX(4)=(2*margProbs(1)-1);% (2*Xsamp(1)-1)
ExpectedValSofX(5)=(2*margProbs(end)-1);% (2*Xsamp(end)-1)


