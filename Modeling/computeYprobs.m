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
%%%%%%%%%%%             Last Modified: 06/15/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the probability distribution of Y=(1/N)*\sum_n X_n
%
% [yProbs,yVals,NMEy] = computeYprobs(p1,transProbs)
%
% Inputs:
%
% p1        
%           Scalar whose value is Pr[X_1=1]
%
% transProbs 
%           (N-1)-by-2 Matrix an such that: 
%           transProbs(n,1) is Pr[x_{n+1}=0|x_{n}=0]
%           transProbs(n,2) is Pr[x_{n+1}=0|x_{n}=1]
%
% Outputs:
%
% yProbs    an (N+1)x1 vector with probabilities of Y=0,1/N,...,1
%
% yVals     an (N+1)x1 vector with values of Y 
%
% NMEy      a scalar representing the normalized entropy of the yProbs
%           distribution

function [yProbs,yVals,NMEy] = computeYprobs(p1,transProbs)


%
% Initialize
%

rng('default'); % use common random numbers (same seed will be used everytime)
                % useful in ESI for reducing variance; see Spall
                % "Introduction to Stochastic Search and Optimization" for
                % more details on CRNs. 
                % can change to rng('shuffle') if you want different values
                % but you increase the variance of the ESI

maxMom              = 4;   % number of moments used in MaxEnt Computations
totMCsamps          = 2^17;% number of Monte Carlo samples used in estimation
threshForAnalytical = log2(totMCsamps)+1;
N                   = size(transProbs,1)+1; % number of CpG sites 

transProbs(transProbs<0)=0;transProbs(transProbs>1)=1; % fix numerical errors

%
% Compute Yvals
%
yVals = 0:(1/N):1;

%
% Compute Yprobs
%
if N<=threshForAnalytical % compute analytically
    
    yProbs = zeros(N+1,1);
    
    for xPatNum=0:((2^N)-1)
        %
        % compute probability of this x pattern
        %
        xBits = dec2bin(xPatNum,N); % convert to string with N bit representation of xPatNum
        
        sumOfX = 0;
        
        if (xBits(1)-'0')==1 % if first bit is a 1
            xProb   = p1;
            prevBit = 1;
            sumOfX  = sumOfX+1;
        elseif (xBits(1)-'0')==0
            xProb   = 1-p1;
            prevBit = 0;
        else 
            disp('Error in computing Y probs')
        end
        
        for n=2:N
            %
            % Update probability to transition to next bit
            %
            currBit = xBits(n)-'0';
            
            if (currBit==0)&&(prevBit==0)
                xProb  = xProb*transProbs(n-1,1);    % Pr[x_{n}=0|x_{n-1}=0] Pr[x_{n-1}=0,..x_{1}]
            elseif (currBit==0)&&(prevBit==1)
                xProb  = xProb*transProbs(n-1,2);    % Pr[x_{n}=0|x_{n-1}=1] Pr[x_{n-1}=1,..x_{1}]
            elseif (currBit==1)&&(prevBit==0)
                xProb  = xProb*(1-transProbs(n-1,1));% Pr[x_{n}=1|x_{n-1}=0] Pr[x_{n-1}=0,..x_{1}]
                sumOfX = sumOfX+1;    
            elseif (currBit==1)&&(prevBit==1)
                xProb  = xProb*(1-transProbs(n-1,2));% Pr[x_{n}=1|x_{n-1}=1] Pr[x_{n-1}=1,..x_{1}]
                sumOfX = sumOfX+1;
            else
                disp('Error in computing Y probs')
            end
            
            % update previous bit for next step in loop
            prevBit=currBit;
        end % end loop over N bits
        
        yProbs(sumOfX+1) = yProbs(sumOfX+1) + xProb;
        
    end % end loop over patterns of X
    
    % correct for numerical problems:
    yProbs(yProbs<0)=0;
    yProbs = yProbs./sum(yProbs);

else % N>16, so compute via MaxEnt and Monte Carlo 
    %
    % Calculate moments using M monte carlo samples
    %
    moments = zeros(maxMom,1);
    for m=1:totMCsamps
        
        Xsamp    = exactSampling(p1,transProbs);
        fracMeth = mean(Xsamp);

        for momNum=1:maxMom
            moments(momNum)=moments(momNum)+(fracMeth^momNum);
        end

    end% end monte carlo sample loop

    moments = moments./totMCsamps;

    %
    % Compute the maxEnt model 
    %
    
    [~,p,~] = maxent(moments,yVals);
    p(p<0)  = 0;
    yProbs  = p./sum(p);
    
end % end decision to compute analytically or via MaxEnt
    


%
% NME of Y
%
yProbsNonZero = yProbs(yProbs>0);
entropy = dot(yProbsNonZero,-log2(yProbsNonZero));
                                
NMEy = entropy/log2(N+1);
if NMEy<0 
    NMEy=0;
elseif NMEy>1
    NMEy=1;
end

end % end of function
