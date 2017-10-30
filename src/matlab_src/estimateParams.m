% informME: An information-theoretic pipeline for WGBS data
% Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu), 
% and Jordi Abante (jabante1@jhu.edu)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
% or see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                   EstimateParams.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 09/08/2017                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a list of BAM files (which correspond to the same 
% phenotype) and estimates the parameters of the 1D Ising model that 
% best fits the methylation data associated with a specific genomic 
% region.
%
% USAGE:
%
% regionStruct = EstimateParams(locationPathName,phenoName,...
%                                 DistInRegion,densityInRegion,dataMat)
% 
% INPUTS:
%
% locationPathName
%               A string that specifies the relative path to the genomic 
%               region being modeled. For example: 'chr1/bp1-3000'. 
%
% phenoName
%               A string that specifies the name of the modeled phenotype.
%
% DistInRegion
%               An Nx1 vector of distances to the next CpG site for each 
%               of the N CpG sites in the genomic region being modeled.
%
% densityInRegion
%               An Nx1 vector of densities for each of the N CpG sites in 
%               the genomic region being modeled.
%
% dataMat
%               A dxN matrix of methylation data within the genomic region 
%               being modeled. Each row of this matrix is a single 
%               observation with elements taking values -1,0,1 
%               corresponding to the methylation status of a CpG site 
%               not being observed (-1), being unmethylated (0), or being 
%               methylated (1).
%

function regionStruct = EstimateParams(locationPathName,phenoName,...
                                 DistInRegion,densityInRegion,dataMat)

regionStruct=[];

% The MCS algorithm does not accept parameters to 'objFnToMinimize' so 
% must make these parameters global variables.

global CpGstart CpGend density Dist newMatrixVec

Dist    = double(DistInRegion(1:(end-1))); 
density = densityInRegion;
N       = length(density);

% Check for sufficient data coverage and writeout coverage to data
% location. 

percentCovered = sum(sum(dataMat>-1,1)>0)/N;
depthOfCov     = sum(sum(dataMat>-1))/N;
if (depthOfCov<2.5) || (percentCovered<(2/3)) 
    return; % Insufficient data - do not build statistical model.
end

% Process data matrix to have only contiguous reads. 

try 
    [newMatrixVec, CpGstart, CpGend] = processMatrix(dataMat);
catch ME
    fprintf(2,['Error processing data matrix from' locationPathName '\n']);
    fprintf(2,getReport(ME));
    return;
end

% Perform maximum likelihood estimation using MCS.

try
    [thetabest,fbest,thetamin,fmi,ncloc,flag] = mcs('feval',...
        'objFnToMinimize',[-10;-100;-20;-3;-3],[10;100;20;3;3]);%#ok<ASGLU>
     thetabest = thetabest(:); % Ensure column vector.                                
catch ME
    fprintf(2,['Error in MCS for phenotype ' phenoName ' at location: ' ...
        locationPathName ' with error:\n']);
    fprintf(2,getReport(ME))
    return;
end                                                                       

try
    
    % Solve model.
    
    [An,Cn] = computeAnCn(density,Dist,thetabest); 
    [logZ1,logZ0,logZ] = computeZ(An,Cn); 
    [logZ1tilde,logZ0tilde,~] = computeZtilde(An,Cn); 
    [p0,transProbs] = computeMCtransProbs(An,Cn,logZ1,logZ0,logZ); 

    % Compute 1D marginals.
    
    margProbs = zeros(N,1); %P(X_n=1)
    margProbs(1) = 1-p0;

    for r=2:N
        logMargProb = calcMargProb(int32(r),int32(0),int32(1),logZ1,...
            logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn);
        margProbs(r) = double(exp(logMargProb));
    end

    % Correct numerical errors.
    
    margProbs(margProbs>(1-eps))   = 1-eps;
    margProbs(margProbs<eps)       = eps;
    transProbs(transProbs>(1-eps)) = 1-eps;
    transProbs(transProbs<eps)     = eps;
    
catch
    fprintf(2,['Error in model estimation for phenotype ' phenoName ...
        'at location: ' locationPathName '\n']);
    return;
end

% Store results.

try 
    regionStruct = struct('DistInRegion',DistInRegion,...
                          'densityInRegion',densityInRegion,...
                          'thetabest',thetabest,...
                          'p0',p0,...
                          'transProbs',transProbs,...
                          'margProbs',margProbs,...
                          'logZ',logZ);
catch
    fprintf(2,['Error in creating structure for phenotype ' phenoName ...
        ' at location: ' locationPathName '\n']);
    return;
end
