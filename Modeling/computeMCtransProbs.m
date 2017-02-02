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
% This function computes the non-homogeneous markov chain transition
% probabilities as well as the initial conditions of this chain, which
% together are an alternative represntation of the Ising Model.
%
% [p1,transProbs] = computeMCtransProbs(An,Cnm,logZ1,logZ0,logZ)
%
% INPUTS:
% An is a Nx1 vector of a_n parameters for n=1,...,N
% Cn is a (N-1)x1 vector of C_n parameters for n=1,...,N-1
% logZ1 is a Nx1 vector of values Z_n(x_n) for x_n=1 and n=1,...,N
% logZ0 is a Nx1 vector of values Z_n(x_n) for x_n=0 and n=1,...,N
% logZ is the partition function
%
% OUTPUTS:
% p1 is the marginal probability  P(X_1=0)
% transProbs is a (N-1)x2 matrix
%   -First  column has probability P(X_{n+1}=0|X_n=0) for n=1,2,...,N-1
%   -Second column has probability P(X_{n+1}=0|X_n=1) for n=1,2,...,N-1
%
% Note the joint dist can be calculated from these probabilities:
% P(\bfX=\bfx) = P(X_1)\prod_{n=1}^{N-1} P(X_{n+1|X_n)
%
% More importantly these probabilities can be used to iteratively draw an
% exact sample from the joint distribution without MCMC.

function [p1,transProbs] = computeMCtransProbs(An,Cnm,logZ1,logZ0,logZ)

N = length(An);
d = 200; % number of digits for variable precision arithmetic

transProbs = sym(zeros(N-1,2));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute boundary condition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
p1 = double(vpa(...
                exp( logZ0(1)-logZ )...
                ,d)); % P(X_N=0)

%
% P(x_2=0|x_1=0) = \phi_1(x_1=0,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=0)
%
transProbs(1,1) = vpa(...
                      exp( -An(1)-An(2)+Cnm(1)+logZ0(2)-logZ0(1) )...
                      ,d);

%
% P(x_2=0|x_1=1) = \phi_1(x_1=1,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=1)
%
transProbs(1,2) = vpa(...
                      exp( An(1)-An(2)-Cnm(1)+logZ0(2)-logZ1(1) )...
                      ,d);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use backwards recursion to compute transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for n=2:(N-1) % i = N-1,N-2,...,3,2
    %
    % P(x_{n+1}=0|x_n=0) = \phi_n(x_n=0,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=0)
    %
    transProbs(n,1) = vpa(...
                          exp( -An(n+1)+Cnm(n)+logZ0(n+1)-logZ0(n) )...
                          ,d);
    
    %
    % P(x_{n+1}=0|x_{n}=1) = \phi_n(x_n=1,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=1)
    %
    transProbs(n,2) = vpa(...
                          exp(-An(n+1)-Cnm(n)+logZ0(n+1)-logZ1(n))...
                          ,d);
end

transProbs = double(transProbs);