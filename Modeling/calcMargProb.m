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
% Calculates the marginal probability of a contiguous set of CpG sites
% taking specific values.
%
% logMargProb = calcMargProb(r,s,x_r_rPLUSs,logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn)
%
% Inputs:
%
% r
%           Index of first CpG sites to have included in the marginal
%
% s         
%           Number of CpG sites beyond the rth to include in the marginal
%
% x_r_rPLUSs
%           Binary vector of length s+1 of the values of methylation for
%           which the probability is being calculated
%
% logZ1 
%           a Nx1 vector of values Z_n(x_n) for x_n=1 and n=1,...,N
%
% logZ0  
%           a Nx1 vector of values Z_n(x_n) for x_n=0 and n=1,...,N
% 
% logZ 
%           the partition function
%
% An 
%           a Nx1 vector of a_n parameters for n=1,...,N
% Cn 
%           a (N-1)x1 vector of C_n parameters for n=1,...,N-1
%
% Outputs:
%
% logMargProb
%           The log (base e) of the marginal probability in question

function logMargProb = calcMargProb(r,s,x_r_rPLUSs,logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn)


d = 100; %set number of digits for VPA

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set to log(1/Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

logMargProb=-logZ;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add log[\tilde{Z}_{r}(x_r)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if x_r_rPLUSs(1)>0 %X(r)=1
    logMargProb = vpa(logMargProb + logZ1tilde(r),d);
else               %X(r)=0
    logMargProb = vpa(logMargProb + logZ0tilde(r),d);
end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add Z_{r+s}(x_{r+s})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if x_r_rPLUSs(s+1)>0 %X(r+s)=1
    logMargProb = vpa(logMargProb + logZ1(r+s),d);
else                 %X(r+s)=0
    logMargProb = vpa(logMargProb + logZ0(r+s),d);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add by \sum_{n=r}^{r+s-1} log[ \phi(x_n,x_{n+1}) ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if s>0 % otherwise this product is empty
    if r>1 % no boundary \phi term needed
        for iter=1:s
            n=iter+r-1;
            logMargProb=vpa(...
                            logMargProb + (An(n+1)+Cn(n)*(2*x_r_rPLUSs(iter)-1))...
                                          *(2*x_r_rPLUSs(iter+1)-1)...
                            ,d);
        end
    else % must compute boundary \phi term since r=1
        %n=1;
        logMargProb=vpa(...
                        logMargProb +  (2*x_r_rPLUSs(1)-1)*An(1)...
                                      +(2*x_r_rPLUSs(2)-1)*( An(2) ...
                                                            +Cn(1)*(2*x_r_rPLUSs(1)-1) )...
                        ,d);
        for n=2:s
            logMargProb=vpa(...
                            logMargProb + (An(n+1)+Cn(n)*(2*x_r_rPLUSs(n)-1))...
                                         *(2*x_r_rPLUSs(n+1)-1)...
                            ,d);
        end                  
    end
end

