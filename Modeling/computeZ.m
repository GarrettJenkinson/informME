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
%%%%%%%%%%%  Statistical Model for DNA Methylation Patterns    %%%%%%%%%%%%
%%%%%%%%%%%  Code by: Garrett Jenkinson                        %%%%%%%%%%%%
%%%%%%%%%%%             Last Modified: 05/15/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function recurses through the Z computation on the Ising model.
%
% [logZ1,logZ0,logZ] = computeZ(An,Cn)
%
% Inputs:
% An is a Nx1 vector of Omega*alpha_n parameters for n=1,...,N
% Cn is a (N-1)x1 vector of Omega*\beta_{n,n+1} parameters for n=1,...,N-1
%
% Outputs:
% logZ1 is a Nx1 vector of values log(Z_n(x_n)) for x_n=1 and n=1,...,N
% logZ0 is a Nx1 vector of values log(Z_n(x_n)) for x_n=0 and n=1,...,N
% logZ  is the log partition function
function [logZ1,logZ0,logZ] = computeZ(An,Cn)



d = 200; % number of digits of numerical accuracy in variable precision arithmetic;  
         % ~32 is roughly floating point accuracy, default setting

N = length(An);

logZ1 = sym(zeros(N,1));
logZ0 = sym(zeros(N,1));

An = sym(An);
Cn  = sym(Cn);

%
% Calculate first boundary values
%
logZ1(N) = 0;
logZ0(N) = 0;

%
% Recurse through non-boundary values
%

for n=(N-1):-1:2
    %
    % Z_n(x_n) = \sum_{x_{n+1}=0}^1 \phi_n(x_n,x_{n+1}) Z_{n+1}(x_{n+1})
    %
    % This implies:
    %
    % log(Z_n(x_n)) = log( exp( log(\phi_n(x_n,0))+log(Z_{n+1}(0) ) ...
    %                     +exp( log(\phi_n(x_n,1))+log(Z_{n+1}(1) ))                            
    %
    logZ1(n) = vpa(...
                   log( exp(-An(n+1)-Cn(n)+logZ0(n+1)) ...
                      + exp( An(n+1)+Cn(n)+logZ1(n+1)) )...
                   ,d);
    logZ0(n) = vpa(...
                   log( exp(-An(n+1)+Cn(n)+logZ0(n+1)) ...
                      + exp( An(n+1)-Cn(n)+logZ1(n+1)) )...
                  ,d);
end

%
% Calculate last boundary values
%

logZ1(1) = vpa(...
               log( exp( An(1)-An(2)-Cn(1)+logZ0(2)) ...
                  + exp( An(1)+An(2)+Cn(1)+logZ1(2)) )...
               ,d);
logZ0(1) = vpa(...
               log( exp(-An(1)-An(2)+Cn(1)+logZ0(2)) ...
                  + exp(-An(1)+An(2)-Cn(1)+logZ1(2)) )...
              ,d);



%
% Compute Log Partition Function = log(Z_1(0)+Z_1(1))
%
logZ = vpa(...
           log( exp(logZ0(1)) + exp(logZ1(1)) )...
           ,d);

