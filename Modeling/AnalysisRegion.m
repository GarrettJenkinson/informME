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
%
% This function computes the statistical summaries in a specified genomic 
% region.
%
% Usage:
% regionStruct = AnalysisRegion(localEstStruct,CpGlocs_local,startBP,endBP,subRegionSize)
%
% This function takes as inputs:
%
% localEstStruct
%               A structure that was generated in parameter estimation for
%               this region.
%
% CpGlocs_local
%               A vector of CpG locations within the current region.
%
% startBP
%               The starting base pair index (1-based) along the genome for
%               this region.
%
% endBP
%               The ending base pair index (1-based) along the genome for
%               this region.
%
% subRegionSize
%               A scalar specifying the number of basepairs in the
%               subregion that should be analyzed. Should divide evenly
%               into the length of the region size (implicitly specified by
%               startBP and endBP)
%
% And returns as an output:
%
% regionStruct
%               A structure of results from analyzing the model in this
%               region. See last line of code for details about what all is
%               included in structure.

function regionStruct = AnalysisRegion(localEstStruct,CpGlocs_local,startBP,endBP,subRegionSize)

subRegStartBPlist = startBP:subRegionSize:endBP;

numSubRegions = length(subRegStartBPlist);

boundaries = [0,0.25,0.5,0.75,1];

ProbsFromMC = zeros(numSubRegions,length(boundaries)-1);
isModeled   = zeros(numSubRegions,1);
Nrq         = zeros(numSubRegions,1);
NMEx        = zeros(numSubRegions,1);
NMEy        = zeros(numSubRegions,1);
MedianVals  = zeros(numSubRegions,1);
MeanAn      = zeros(numSubRegions,1);
MeanCn      = zeros(numSubRegions,1);
MeanLevels  = zeros(numSubRegions,1);
NCU         = zeros(numSubRegions,1);
NRDE        = zeros(numSubRegions,1);
ESI         = zeros(numSubRegions,1);
NNcorr      = zeros(numSubRegions,1);
TURN        = zeros(numSubRegions,1);
MargEnt     = zeros(numSubRegions,1);

% sparse matrix parameters
rows=nan(floor(subRegionSize/2)*numSubRegions,1);
cols=nan(floor(subRegionSize/2)*numSubRegions,1);
vals=nan(floor(subRegionSize/2)*numSubRegions,1);
sparseIndex = 1;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize model from data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
normEntReg      = localEstStruct.normEntRegion;  
thetabest       = localEstStruct.thetabest;
margProbs       = localEstStruct.margProbs;
DistInRegion    = double(localEstStruct.DistInRegion);
densityInRegion = localEstStruct.densityInRegion;
transProbs      = localEstStruct.transProbs;
NNcorrCpG       = localEstStruct.NNcorr;
[An,Cn]         = computeAnCnm(densityInRegion,DistInRegion(1:(end-1)),thetabest);

%correct numerical errors:
margProbs(margProbs<0)=0;
margProbs(margProbs>1)=1;

% Compute marginal entropies
margEntCpGs = h_func(margProbs);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve "differential" models for sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
epsilon=0.01;
N=length(An);

[Ana,Cna]         = computeAnCnm(densityInRegion,DistInRegion(1:(end-1)),thetabest*(1+epsilon));

[logZ1a,logZ0a,logZa] = computeZ(Ana,Cna); 
logZ1a=vpa(logZ1a);logZ0a=vpa(logZ0a);logZa =vpa(logZa); % make output syms
[logZ1tildea,logZ0tildea,~] = computeZtilde(Ana,Cna); 
logZ1tildea=vpa(logZ1tildea);logZ0tildea=vpa(logZ0tildea); %make output syms
[p0a,transProbsa] = computeMCtransProbs(Ana,Cna,logZ1a,logZ0a,logZa); 



% Compute 1D marginals
margProbsa = zeros(N,1); %P(X_n=1)
margProbsa(1) = 1-p0a;


for r=2:N
    %s = 0; %x_r_rPLUSs = 1;
    logMargProba = calcMargProb(int32(r),int32(0),int32(1),logZ1a,logZ0a,logZa,...
                               logZ1tildea,logZ0tildea,Ana,Cna);
    margProbsa(r) = double(exp(logMargProba));
end

%
% correct numerical errors
%
margProbsa(margProbsa>(1-eps))   = 1-eps;
margProbsa(margProbsa<eps)       = eps;
transProbsa(transProbsa>(1-eps)) = 1-eps;
transProbsa(transProbsa<eps)     = eps;


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

subRegCount=0;

for subRegStartBP = subRegStartBPlist
    
    subRegCount = subRegCount+1;
    subRegEndBP = int64(subRegStartBP+subRegionSize-1);
    
    % find CpG sites 
    [lower_index,upper_index] = findSortedIndices(CpGlocs_local,subRegStartBP,subRegEndBP);
    CpGlocs                   = CpGlocs_local(lower_index:upper_index);
    Nrq(subRegCount)          = length(CpGlocs);
    
    
    % calc probs via Monte Carlo in subregion
    if Nrq(subRegCount)>0
        
        isModeled(subRegCount) = 1; % set flag, this region is modeled
        
        %
        % Compute mean values
        %
        MeanAn(subRegCount)     = mean(An(lower_index:upper_index));
        MeanCn(subRegCount)     = mean(Cn(max(lower_index-1,1):min(upper_index,length(Cn))));
        NNcorr(subRegCount)     = mean(NNcorrCpG(max(lower_index-1,1):min(upper_index,length(NNcorrCpG))));
        MeanLevels(subRegCount) = mean(margProbs(lower_index:upper_index));
        MargEnt(subRegCount)    = mean(margEntCpGs(lower_index:upper_index));
        
        %
        % Compute Capacities and TURN variables
        %
       
        lambdaVals  = margProbs(lower_index:upper_index)./(1-margProbs(lower_index:upper_index)); % P(1)/P(0)
        
        TURN(subRegCount) = mean(log2(lambdaVals));
        

        % Compute Capacity
        NCUvals = zeros(size(lambdaVals));
        
        NCUvals(lambdaVals>=1) = 1 - 0.52.*h_func( lambdaVals(lambdaVals>=1)./(1+lambdaVals(lambdaVals>=1)) )...
                                        ./(1+lambdaVals(lambdaVals>=1));
        NCUvals(lambdaVals<1) = 1 - 0.52.*h_func( lambdaVals(lambdaVals<1)./(1+lambdaVals(lambdaVals<1)) )...
                                     .*lambdaVals(lambdaVals<1)./(1+lambdaVals(lambdaVals<1));

        NCU(subRegCount) = mean(NCUvals);

        % Compute NRDE
        NRDEvals = zeros(size(lambdaVals));

        NRDEvals(lambdaVals<1) = log2( (1+lambdaVals(lambdaVals<1))...
                                      ./(2*lambdaVals(lambdaVals<1)) ) ...
                                 + 4.76;
        NRDEvals(lambdaVals>=1) = log2( (1+lambdaVals(lambdaVals>=1))./2 ) ...
                                  + 4.76;

        NRDE(subRegCount) = mean(NRDEvals);

        
        %
        % Compute values based on Y=\sum_n X_n
        %
        if Nrq(subRegCount) >= 2 % 2 or more CpGs
            
            [pNorm,yVal,NMEyTemp] = computeYprobs(margProbs(lower_index),transProbs(lower_index:(upper_index-1),:));
            
            if sum(isnan(pNorm))>0
                isModeled(subRegCount) = 0;
                continue;
            end
            
            NMEy(subRegCount)=NMEyTemp;
            
            %
            % Add probability to sparse matrix:
            % fullProbs(subRegCount,1:length(pNorm))=pNorm;
            %
            for indexNum=1:length(pNorm)
                vals(sparseIndex)=pNorm(indexNum);
                rows(sparseIndex)=subRegCount;
                cols(sparseIndex)=indexNum;
                sparseIndex=sparseIndex+1;
            end

            
            %
            % Compute coarse Y probabilities
            %
            
            if ~isempty(pNorm(yVal==boundaries(3)))
                correction = (pNorm(yVal==boundaries(3))/2);
            else
                correction = 0;
            end
            
            ProbsFromMC(subRegCount,1) = sum(pNorm((boundaries(1)<=yVal)&(yVal<=boundaries(2))));% 0<= Y <= 0.25
            ProbsFromMC(subRegCount,2) = sum(pNorm((boundaries(2)<yVal)&(yVal<boundaries(3)))) + correction; % 0.25< Y <= 0.5
            ProbsFromMC(subRegCount,3) = sum(pNorm((boundaries(3)<yVal)&(yVal<boundaries(4)))) + correction;  % 0.5<= Y < 0.75
            ProbsFromMC(subRegCount,4) = sum(pNorm((boundaries(4)<=yVal)&(yVal<=boundaries(5)))); % 0.75=< Y <= 1
            
            
            %
            % Find medians
            %
            med_index = find(cumsum(pNorm)>=0.5,1,'first');
            strengthVals = 0:Nrq(subRegCount);
            MedianVals(subRegCount) = strengthVals(med_index);
            
            %
            % Do sensitivity calculation
            %
            [~,~,NMEa] = computeYprobs(margProbsa(lower_index),transProbsa(lower_index:(upper_index-1),:));
                        
            if sum(isnan(NMEa))==0 
                %
                % Compute JS distance
                %

%                 pNorm  = pNorm(:);  % ensures column vector
%                 pNorma = pNorma(:); % ensures column vector
%                 

%                 % JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
%                 % where M=(P+M)/2
%                 % and KL(P|M)=sum_i P_i log2(P_i/M_i)
%                 Dpqa = sqrt( sum(   pNorm(pNorm>0).*log2(2.*pNorm(pNorm>0)  ./(pNorm(pNorm>0) +pNorma(pNorm>0)))  )/2 ...
%                             +sum( pNorma(pNorma>0).*log2(2.*pNorma(pNorma>0)./(pNorm(pNorma>0)+pNorma(pNorma>0))) )/2 );
                Dpqa=abs(NMEa-NMEy(subRegCount))*log2(Nrq(subRegCount)+1);

                ESI(subRegCount) = Dpqa/epsilon;
                
            end
              
        else% Nrq(subRegCount)==1
            %
            % just use exact marginal probabilities
            %
            ProbsFromMC(subRegCount,1)=1-margProbs(lower_index);
            ProbsFromMC(subRegCount,4)=margProbs(lower_index);
            
            NMEx(subRegCount)= -margProbs(lower_index).*log2(margProbs(lower_index))...
                           -(1-margProbs(lower_index)).*log2(1-margProbs(lower_index));
             
            NMEy(subRegCount) = NMEx(subRegCount); 
            
            MedianVals(subRegCount,1) = margProbs(lower_index)>0.5; % zero if P(X=0)>=0.5, one else
            
            %
            % Add probs to sparse matrix:
            % fullProbs(subRegCount,1:2)=[1-margProbs(lower_index),margProbs(lower_index)];
            %
            vals(sparseIndex)=1-margProbs(lower_index);
            vals(sparseIndex+1)=margProbs(lower_index);
            rows(sparseIndex)=subRegCount;
            rows(sparseIndex+1)=subRegCount;
            cols(sparseIndex)=1;
            cols(sparseIndex+1)=2;
            sparseIndex=sparseIndex+2;
            
            %
            % Do sensitivity calculation
            %
            
%            pNorm  = [1-margProbs(lower_index);margProbs(lower_index)];  
            pNorma = [1-margProbsa(lower_index);margProbsa(lower_index)];  
            NMEa = -pNorma'*log2(pNorma);
            
%             % JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
%             % where M=(P+M)/2
%             % and KL(P|M)=sum_i P_i log2(P_i/M_i)
%             Dpqa = sqrt( sum(   pNorm(pNorm>0).*log2(2.*pNorm(pNorm>0)  ./(pNorm(pNorm>0) +pNorma(pNorm>0)))  )/2 ...
%                         +sum( pNorma(pNorma>0).*log2(2.*pNorma(pNorma>0)./(pNorm(pNorma>0)+pNorma(pNorma>0))) )/2 );
            
            Dpqa = abs(NMEa-NMEy(subRegCount));
            
            ESI(subRegCount) = Dpqa/epsilon;
                
        end % end  if Nrq(subRegCount) >=2 else Nrq(subRegCount)==1

    end% end Nrq(subRegCount)>0 condition
    
end % end loop over subregions
 
%
% Create sparse matrix fullProbs with rows containing the probability
% distribution of Y
%
sparseIndex=sparseIndex-1;
fullProbs = sparse(rows(1:sparseIndex),cols(1:sparseIndex),vals(1:sparseIndex),numSubRegions,max(cols(1:sparseIndex)));

%
% Save all results into a structure
%
regionStruct = struct('ProbsFromMC',double(ProbsFromMC),'isModeled',double(isModeled),'Nrq',double(Nrq),...
                       'NMEy',double(NMEy),'NMEx',double(NMEx),'MedianVals',double(MedianVals),...
                       'fullProbs',fullProbs,'MeanAn',MeanAn,'MeanCn',MeanCn,'MeanLevels',MeanLevels,...
                       'CpGlocs_local',CpGlocs_local,'ESI',ESI,'NCU',NCU,'NRDE',NRDE,'thetabest',thetabest,...
                       'NNcorr',NNcorr,'TURN',TURN,'normEntReg',normEntReg,'MargEnt',MargEnt);

end
