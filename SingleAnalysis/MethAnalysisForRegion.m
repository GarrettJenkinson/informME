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
%%%%%%%%                 MethAnalysisForRegion.m                   %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/02/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs methylation analysis of a genomic region used 
% for estimating the parameters of the Ising model by computing a number 
% of statistical summaries of the methylation state within the region, 
% including probability distributions of methylation levels, mean 
% meathylation levels, and normalized methylation entropies. If desired, 
% this function also computes entropic sensitivity indices, as well 
% information-theoretic quantities associated with methylation channels, 
% such as trunover ratios, channel capacities, and relative dissipated 
% energies. 
%
% USAGE:
%
% regionStruct = MethAnalysisForRegion(localEstStruct,CpGlocs_local,...
%                              startBP,endBP,subregionSize,ESIflag,MCflag)
%
% INPUTS:
%
% localEstStruct
%               A structure generated during parameter estimation within 
%               the genomic region.
%
% CpGlocs_local
%               A vector of CpG locations within the genomic region.
%
% startBP
%               The starting base pair index (1-based) along the genome for
%               the genomic region.
%
% endBP
%               The ending base pair index (1-based) along the genome for
%               the genomic region.
%
% subregionSize
%               A scalar that specifies the number of base pairs within 
%               subregions of the genomic region determining the 
%               resolution of methylation analysis. The analysis subregions 
%               must be of the same length (in base pairs) and chosen 
%               so that the genomic region is partitioned into an integer 
%               number of nonoverlapping subregions.
%
% ESIflag
%               Flag that determines whether this function performs 
%               computation of the entropic sensitivity index (ESI). 
%               0: no ESI computation. 
%               1: allow ESI computation.
%               Default value: 0
%
% MCflag
%               Flag that determines whether this function performs 
%               computation of turnover ratios, CpG entropies, capacities, 
%               and relative dissipated energies of methylation 
%               channels (MCs). 
%               0: no MC computations. 
%               1: allow MC computations. 
%               Default value: 0
%
% OUTPUT:
%
% regionStruct
%               A structure summarizing the methylation analysis results 
%               containing the following information:
%               o The locations of the CpG sites within the genomic region.
%               o Numbers of CpG sites within the analysis subregions. 
%               o Which analysis subregions are modeled and which are not.
%               o Estimated parameters of Ising model in genomic region
%               o Methylation level probabilities in modeled subregions.
%               o Coarse methylation level probabilities.
%               o Mean methylation levels.
%               o Normalized methylation entropies.
%               o Entropic sensitivity indices (if ESIflag = 1).
%               o Turnover ratios (if MCflag = 1).
%               o Channel capacities (if MCflag = 1).
%               o Relative dissipated energies (if MCflag = 1).
%

function regionStruct = MethAnalysisForRegion(localEstStruct,CpGlocs_local,...
							   startBP,endBP,subregionSize,ESIflag,MCflag)

% Initialization.

subRegStartBPlist = startBP:subregionSize:endBP; 
                              % Start base pairs for analysis subregions.
numSubRegions     = length(subRegStartBPlist); 
                              % Number of analysis subregions.
isModeled         = zeros(numSubRegions,1); 
                              % Vector of 1's and 0's indicating whether a 
                              % subregion is modeled (1) or not modeled (0).
thresh            = [0,0.25,0.5,0.75,1]; 
                              % Thresholds for determing coarse methylation 
                              % level probabilities.
cLProbs           = zeros(numSubRegions,length(thresh)-1); 
                              % Coarse methylation level probabilities.
Ncg               = zeros(numSubRegions,1); 
                              % Number of CpGs in each analysis subregion.
MML               = zeros(numSubRegions,1); 
                              % Mean methylation level.
NME               = zeros(numSubRegions,1); 
                              % Normalized methylation entropy.
                     
if ESIflag
	ESI = zeros(numSubRegions,1); % Entropic sensitivity index.
else
	ESI = [];
end

if MCflag
    TURN = zeros(numSubRegions,1); % Turnover ratio.
    CAP  = zeros(numSubRegions,1); % Channel capacity.
	RDE  = zeros(numSubRegions,1); % Relative dissipated energy.
else
    TURN = [];
	CAP  = [];
	RDE  = [];

end

% Sparse matrix parameters. Used to store methylation level 
% probabilities in analysis subregions.
rows = nan(floor(subregionSize/2)*numSubRegions,1);
cols = nan(floor(subregionSize/2)*numSubRegions,1);
vals = nan(floor(subregionSize/2)*numSubRegions,1);
sparseIndex = 1;

% Model from data structure.
thetabest       = localEstStruct.thetabest;
margProbs       = localEstStruct.margProbs;
DistInRegion    = double(localEstStruct.DistInRegion);
densityInRegion = localEstStruct.densityInRegion;
transProbs      = localEstStruct.transProbs;

% Correct numerical errors.
margProbs(margProbs<0) = 0;
margProbs(margProbs>1) = 1;

if ESIflag % Solve "differential" models for sensitivity analysis.
	
	epsilon = 0.01;
	N = length(CpGlocs_local);

	[Ana,Cna] = computeAnCn(densityInRegion,DistInRegion(1:(end-1)), ...
							thetabest*(1+epsilon));

	[logZ1a,logZ0a,logZa] = computeZ(Ana,Cna); 
	logZ1a = vpa(logZ1a);
    logZ0a = vpa(logZ0a);
    logZa  = vpa(logZa);
	[logZ1tildea,logZ0tildea,~] = computeZtilde(Ana,Cna); 
	logZ1tildea = vpa(logZ1tildea);
    logZ0tildea = vpa(logZ0tildea);
	[p0a,transProbsa] = computeMCtransProbs(Ana,Cna,logZ1a,logZ0a,logZa); 

	% Compute 1D marginal probabilities.
    
	margProbsa    = zeros(N,1); % Pr[X_n=1]
	margProbsa(1) = 1-p0a;

	for r = 2:N
		logMargProba = calcMargProb(int32(r),int32(0),int32(1),logZ1a,...
                             logZ0a,logZa,logZ1tildea,logZ0tildea,Ana,Cna);
		margProbsa(r) = double(exp(logMargProba));
	end

	% Correct numerical errors,
    
	margProbsa(margProbsa>(1-eps))   = 1-eps;
	margProbsa(margProbsa<eps)       = eps;
	transProbsa(transProbsa>(1-eps)) = 1-eps;
	transProbsa(transProbsa<eps)     = eps;
    
end % End ESIflag.

% Loop through analysis subregions.

subRegCount = 0;

for subRegStartBP = subRegStartBPlist
    
    subRegCount = subRegCount+1;
    subRegEndBP = int64(subRegStartBP+subregionSize-1);
    
    % Find CpG sites. 
    
    [lower_index,upper_index] = findSortedIndices(CpGlocs_local,...
											subRegStartBP,subRegEndBP);
    CpGlocs          = CpGlocs_local(lower_index:upper_index);
    Ncg(subRegCount) = length(CpGlocs);
    
    % Calculate methylation level probabilities in subregion. 
    
    if Ncg(subRegCount) > 0
        
        isModeled(subRegCount) = 1; % Set flag - this subregion is modeled.
        
        % Compute mean methylation level within subregion.
        
        MML(subRegCount) = mean(margProbs(lower_index:upper_index));
       
        % Compute normalized methylation entropy.

        if Ncg(subRegCount) >= 2 % Two or more CpGs.
           
            [LProbs,LVals,NMETemp] = computeLstats(margProbs(lower_index),...
								transProbs(lower_index:(upper_index-1),:));
            
            if sum(isnan(LProbs)) > 0
                isModeled(subRegCount) = 0;
                continue;
            end
            
            NME(subRegCount) = NMETemp;
            
            % Add probabilities to sparse matrix 
            % fullLProbs(subRegCount,1:length(LProbs)) = LProbs. 
      
            for indexNum=1:length(LProbs)
                vals(sparseIndex) = LProbs(indexNum);
                rows(sparseIndex) = subRegCount;
                cols(sparseIndex) = indexNum;
                sparseIndex       = sparseIndex+1;
            end

            % Compute coarse methylation level probabilities. 
            
            if ~isempty(LProbs(LVals==thresh(3)))
                correction = (LProbs(LVals==thresh(3))/2);
            else
                correction = 0;
            end

            % Pr[0 <= L <= 0.25]
            cLProbs(subRegCount,1) = sum(LProbs((thresh(1)<=LVals)&...
									       (LVals<=thresh(2))));
			% Pr[0.25 < L < 0.5] + 0.5*Pr[L = 0.5]
            cLProbs(subRegCount,2) = sum(LProbs((thresh(2)<LVals)&...
									       (LVals<thresh(3)))) + correction; 
			% 0.5*Pr[L = 0.5] + Pr[0.5 < L < 0.75]
            cLProbs(subRegCount,3) = sum(LProbs((thresh(3)<LVals)&...
									       (LVals<thresh(4)))) + correction; 
			% Pr[0.75 <= L <= 1] 
            cLProbs(subRegCount,4) = sum(LProbs((thresh(4)<=LVals)&...
                                           (LVals<=thresh(5))));
                                       
            if ESIflag
                % Do sensitivity calculation.
                
                [~,~,NMEa] = computeLstats(margProbsa(lower_index),...
							   transProbsa(lower_index:(upper_index-1),:));
		                    
                if sum(isnan(NMEa))==0
                    Da = abs(NMEa-NME(subRegCount))*log2(Ncg(subRegCount)+1);
                    ESI(subRegCount) = Da/epsilon;
                end
                
            end
              
        else % Ncg(subRegCount) == 1
            
            % Just use exact marginal probabilities.
            
            cLProbs(subRegCount,1) = 1-margProbs(lower_index);
            cLProbs(subRegCount,4) = margProbs(lower_index);
            
            NME(subRegCount) = ...
                -margProbs(lower_index).*log2(margProbs(lower_index))...
                -(1-margProbs(lower_index)).*log2(1-margProbs(lower_index));
            
            % Add methylation level probabilities to sparse matrix 
            % fullLProbs(subRegCount,1:2) = 
            %            [1-margProbs(lower_index),margProbs(lower_index)];
          
            vals(sparseIndex)   = 1-margProbs(lower_index);
            vals(sparseIndex+1) = margProbs(lower_index);
            rows(sparseIndex)   = subRegCount;
            rows(sparseIndex+1) = subRegCount;
            cols(sparseIndex)   = 1;
            cols(sparseIndex+1) = 2;
            sparseIndex=sparseIndex+2;
            
            if ESIflag
                
                % Do sensitivity calculation.
                
                LProbsa = [1-margProbsa(lower_index);margProbsa(lower_index)];
                NMEa    = -LProbsa'*log2(LProbsa);
                Da      = abs(NMEa-NME(subRegCount));
                ESI(subRegCount) = Da/epsilon;
                
            end
                
        end % End if Ncg(subRegCount) >= 2.
        
        if MCflag
            
            % Compute average turnover ratio within subregion.
            
            lambdaVals = margProbs(lower_index:upper_index)...
                             ./(1-margProbs(lower_index:upper_index));
            TURN(subRegCount) = mean(log2(lambdaVals));
            
            % Compute average capacity within subregion.
            
            CAPvals = zeros(size(lambdaVals));
            CAPvals(lambdaVals>=1) = 1-0.52.*h_func(lambdaVals(lambdaVals>=1)...
              ./(1+lambdaVals(lambdaVals>=1)))./(1+lambdaVals(lambdaVals>=1));
            CAPvals(lambdaVals<1) = 1-0.52.*h_func(lambdaVals(lambdaVals<1)...
                ./(1+lambdaVals(lambdaVals<1))).*lambdaVals(lambdaVals<1)...
                ./(1+lambdaVals(lambdaVals<1));
            CAP(subRegCount) = mean(CAPvals);
            
            % Compute average relative dissipated energy within subregion.
            
            RDEvals = zeros(size(lambdaVals));
            RDEvals(lambdaVals<1) = log2((1+lambdaVals(lambdaVals<1))...
                                    ./(2*lambdaVals(lambdaVals<1))) + 4.76;
            RDEvals(lambdaVals>=1) = ...
                             log2((1+lambdaVals(lambdaVals>=1))./2) + 4.76;
            RDE(subRegCount) = mean(RDEvals);
            
        end % End MCflag

    end % End Ncg(subRegCount) > 0 condition.
    
end % End loop over analysis subregions.
 
% Create sparse matrix fullProbs with rows containing the probability
% distribution of L. 

sparseIndex = sparseIndex-1;
fullLProbs  = sparse(rows(1:sparseIndex),cols(1:sparseIndex),...
                           vals(1:sparseIndex),numSubRegions,...
                           max(cols(1:sparseIndex)));

% Save all results into a structure. 

regionStruct = struct('CpGlocs_local',CpGlocs_local,...
                      'Ncg',double(Ncg),...
                      'isModeled',double(isModeled),...
                      'thetabest',thetabest,...
                      'fullLProbs',fullLProbs,...
                      'cLProbs',double(cLProbs),...
                      'MML',MML,...
                      'NME',double(NME),...
                      'ESI',ESI,...
                      'TURN',TURN,...
					  'CAP',CAP,...
					  'RDE',RDE);

end