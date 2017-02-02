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
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                    computeLstats.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/07/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the probability distribution of the methylation 
% level L within a genomic region comprised N CpG sites as well as the 
% corresponding normalized methylation entropy.
%
% USAGE:
%
% [LProbs,LVals,NME] = computeLstats(p1,transProbs)
%
% INPUTS:
%
% p1        
%           The probability Pr[X_1=1] of the first CpG site in the genomic 
%           region to be methylated. 
%
% transProbs 
%           (N-1)x2 matrix whose elements are given by:  
%               transProbs(n,1) = Pr[X_{n+1}=0 | X_n=0]
%               transProbs(n,2) = Pr[X_{n+1}=0 | X_n=1]
%            where X_n is the methylation status of the n-th CpG 
%            site within the genomic region containing N CpG sites.
%
% OUTPUTS:
%
% LProbs    
%           (N+1)x1 vector containing the probabilities of the 
%            methylation level L to take values 0,1/N,...,1.
%
% LVals     
%           (N+1)x1 vector containing the methylation level 
%            values 0,1/N,...,1.
%
% NME       The normalized methylation entropy.
%

function [LProbs,LVals,NME] = computeLstats(p1,transProbs)

% Initialization.

rng('default'); % Use common random numbers (same seed will be used each 
                % time). Useful for reducing estimation variance. Can be  
                % changed to rng('shuffle') if different values are
                % desirable, but this increases variance.

maxMom              = 4; % Number of moments used in MaxEnt computations.
totMCsamps          = 2^17; % Number of Monte Carlo samples used in estimation.
threshForAnalytical = log2(totMCsamps)+1; % Threshold for analytical computations.
N                   = size(transProbs,1)+1; % Number of CpG sites.  

transProbs(transProbs<0)=0;transProbs(transProbs>1)=1; % Fix numerical errors.

% Compute Lvals.

LVals = 0:(1/N):1;

% Compute Lprobs.

if N <= threshForAnalytical % Compute analytically.
    
    LProbs = zeros(N+1,1);
    
    for xPatNum=0:((2^N)-1)
        
        % Compute probability of this xPatNum pattern.
        
        xBits = dec2bin(xPatNum,N); % Convert to string with N bit 
                                    % representation of xPatNum.
        
        sumOfX = 0;
        
        if (xBits(1)-'0')==1 % If first bit is 1.
            xProb   = p1;
            prevBit = 1;
            sumOfX  = sumOfX+1;
        elseif (xBits(1)-'0')==0
            xProb   = 1-p1;
            prevBit = 0;
        else 
            disp('Error in computing methylation level probabilities')
        end
        
        for n = 2:N
            
            % Update probability of transition to next bit.
            
            currBit = xBits(n)-'0';
            
            if (currBit==0)&&(prevBit==0)
                xProb = xProb*transProbs(n-1,1); 
            elseif (currBit==0)&&(prevBit==1)
                xProb  = xProb*transProbs(n-1,2);
            elseif (currBit==1)&&(prevBit==0)
                xProb  = xProb*(1-transProbs(n-1,1));
                sumOfX = sumOfX+1;    
            elseif (currBit==1)&&(prevBit==1)
                xProb  = xProb*(1-transProbs(n-1,2));
                sumOfX = sumOfX+1;
            else
                disp('Error in computing methylation level probabilities')
            end
            
            prevBit=currBit; % Update previous bit for next step in loop.
        end % End loop over N bits.
        
        LProbs(sumOfX+1) = LProbs(sumOfX+1) + xProb;
        
    end % End loop over xPatNum patterns. 
    
    % Correct for numerical problems.
    
    LProbs(LProbs<0) = 0;
    LProbs = LProbs./sum(LProbs);

else
    
    % Compute via Monte Carlo and maxEnt.
    
    % Estimate moments using Monte Carlo.
    
    moments = zeros(maxMom,1);
    
    for m = 1:totMCsamps
        
        Xsamp    = exactSampling(p1,transProbs);
        fracMeth = mean(Xsamp);

        for momNum=1:maxMom
            moments(momNum)=moments(momNum)+(fracMeth^momNum);
        end

    end % End Monte Carlo loop.

    moments = moments./totMCsamps;

    % Compute maxEnt model. 
  
    [~,p,~] = maxent(moments,LVals);
    p(p<0)  = 0;
    LProbs  = p./sum(p);
    
end % End decision to compute analytically or via maxEnt.

% Compute normalized methylation entropy.

LProbsNonZero = LProbs(LProbs>0);
entropy = dot(LProbsNonZero,-log2(LProbsNonZero));
                                
NME = entropy/log2(N+1);
if NME < 0 
    NME = 0;
elseif NME>1
    NME = 1;
end

end
