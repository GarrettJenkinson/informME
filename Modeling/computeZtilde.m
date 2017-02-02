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
% This function recurses through the Ztilde computation on the Ising model
%
% [logZ1tilde,logZ0tilde,logZtilde] = computeZtilde(An,Cn)
%
% Inputs:
%
% An is a Nx1 vector of Omega*alpha_n parameters for n=1,...,N
%
% Cn is a (N-1)x1 vector of Omega*\beta_{n,n+1} parameters for n=1,...,N-1
%
% Outputs:
%
% logZ1tilde is a Nx1 vector of values log(\tilde{Z}_n(x_n)) for x_n=1 and n=1,...,N
%
% logZ0tilde is a Nx1 vector of values log(\tilde{Z}_n(x_n)) for x_n=0 and n=1,...,N
%
% logZtilde is the log partition function

function [logZ1tilde,logZ0tilde,logZtilde] = computeZtilde(An,Cn)



N = length(An);

d = 100; % number of digits of numerical accuracy in variable precision arithmetic;  
        % ~32 is roughly floating point accuracy, default setting

logZ1tilde = sym(zeros(N,1));
logZ0tilde = sym(zeros(N,1));


An   = sym(An);
Cn  = sym(Cn);


%
% Calculate first two boundary values
%
logZ1tilde(1) = 0;
logZ0tilde(1) = 0;
logZ1tilde(2) = vpa(...
                    log( exp(-An(1)+An(2)-Cn(1) ) ...
                        +exp( An(1)+An(2)+Cn(1) ) )...
                    ,d);
logZ0tilde(2) = vpa(...
                    log( exp(-An(1)-An(2)+Cn(1) ) ...
                        +exp( An(1)-An(2)-Cn(1) ) )...
                    ,d);

%
% Recurse through non-boundary values
%
for n=3:N
    %
    % \tilde{Z}_n(x_n) = \sum_{x_{n-1}=0}^1 \phi_{n-1}(x_{n-1},x_n) \tilde{Z}_{n-1}(x_{n-1})
    %
    % This implies:
    %
    % log(\tilde{Z}_n(x_n)) = log( exp( log(\phi_{n-1}(0,x_n))+log(Z_{n-1}(0) ) ...
    %                             +exp( log(\phi_{n-1}(1,x_n))+log(Z_{n-1}(1) ))                            
    %
    logZ1tilde(n) = vpa(...
                        log( exp( An(n)-Cn(n-1)+logZ0tilde(n-1) ) ...
                            +exp( An(n)+Cn(n-1)+logZ1tilde(n-1) ) )...
                        ,d);
    logZ0tilde(n) = vpa(...
                        log( exp(-An(n)+Cn(n-1)+logZ0tilde(n-1) ) ...
                            +exp(-An(n)-Cn(n-1)+logZ1tilde(n-1) ) )...
                        ,d);
end

%
% Compute Log Partition Function = log(Z_N(0)+Z_N(1))
%
logZtilde = vpa(...
                log( exp(logZ0tilde(N)) + exp(logZ1tilde(N)) )...
                ,d);

