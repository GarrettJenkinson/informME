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
%
% This function takes a list of bamFileNames that come from the same
% phenotype and fits the statistical model for methylation patterns in the
% specified genomic region.
%
% regionStruct = EstimateParams(locationPathName,phenoName,DistInRegion,densityInRegion,dataMat)
%
% This function takes as inputs:
%
% locationPathName
%               A string with the relative path to the location being
%               modeled. For example: 'chr1/bp1-3000' 
%
% phenoName
%               A string with the name of the phenotype with data stored in
%               the BAM files specified in bamFileNames.
%
% DistInRegion
%               A Nx1 vector of distances to the next CpG site for each CpG
%               site in the region
%
% densityInRegion
%               A Nx1 vector of densities for each CpG site in the region
%
% dataMat
%               A dxN matrix of data. Each row is a single observation with
%               elements taking values -1,0,1 corresponding to the CpG site
%               being unobserved, being unmethylated, or being methylated, 
%               respectively.

function regionStruct = EstimateParams(locationPathName,phenoName,DistInRegion,densityInRegion,dataMat)


regionStruct=[];

% MCS algorithm does not accept parameters to 'objFnToMinimize
% so must make these parameters global variables
global CpGstart CpGend density Dist newMatrix 


Dist             = double(DistInRegion(1:(end-1))); 
density          = densityInRegion;
N                = length(density);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for sufficient data coverage and writeout coverage to data location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

percentCovered = sum(sum(dataMat>-1,1)>0)/N;
depthOfCov     = sum(sum(dataMat>-1))/N;
if (depthOfCov<2.5) || (percentCovered<(2/3)) 
    return; % insufficient data; do not build statistical model
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process data matrix to have only contiguous reads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
try 
    [newMatrix, CpGstart, CpGend] = processMatrix( dataMat );
catch
    disp(['Error processing data matrix from: ' dataFile]);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute maximum likelihood estimator via MCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

try
    [thetabest,fbest,thetamin,fmi,ncloc,flag]=mcs('feval','objFnToMinimize',...
                                          [-10;-100;-20;-3;-3],[10;100;20;3;3]); %#ok<ASGLU>
     thetabest=thetabest(:); %ensure column vector                                
catch ME
    display(['Error in MCS for phenotype ' phenoName ' at location: ' locationPathName ' with stack:']);
    getReport(ME)
    ME.stack.file
    ME.stack.name
    ME.stack.line
    ME.identifier
    return;
end                                                                       


try
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %

    [invOmega_An,invOmega_Cnm] = computeAnCnm(density,Dist,thetabest); 
    [logZ1,logZ0,logZ] = computeZ(invOmega_An,invOmega_Cnm); 
    logZ1=vpa(logZ1);logZ0=vpa(logZ0);logZ =vpa(logZ); % make output syms
    [logZ1tilde,logZ0tilde,~] = computeZtilde(invOmega_An,invOmega_Cnm); 
    logZ1tilde=vpa(logZ1tilde);logZ0tilde=vpa(logZ0tilde); %make output syms
    [p0,transProbs] = computeMCtransProbs(invOmega_An,invOmega_Cnm,logZ1,logZ0,logZ); 

    %
    % Compute 1D marginals
    %
    margProbs = zeros(N,1); %P(X_n=1)
    margProbs(1) = 1-p0;

    for r=2:N
        %s = 0; %x_r_rPLUSs = 1;
        logMargProb = calcMargProb(int32(r),int32(0),int32(1),logZ1,logZ0,logZ,...
                                   logZ1tilde,logZ0tilde,invOmega_An,invOmega_Cnm);
        margProbs(r) = double(exp(logMargProb));
    end

    %
    % correct numerical errors
    %
    margProbs(margProbs>(1-eps))   = 1-eps;
    margProbs(margProbs<eps)       = eps;
    transProbs(transProbs>(1-eps)) = 1-eps;
    transProbs(transProbs<eps)     = eps;

    %
    % Compute ExpectedValSofX and entropy for the region
    %
    [ExpectedValSofX,NNcorr] = expectValSXexact(margProbs,transProbs,density,Dist); 
    ExpectedValSofX=ExpectedValSofX(:);%ensure column vector
    
    entrRegion    = logZ - dot(thetabest,ExpectedValSofX);
    normEntRegion = entrRegion/(N*log(2));  
    
    if normEntRegion>1 || normEntRegion<0
        normEntRegion=nan;
    end
    
catch
    disp(['Error in Model Solving for phenotype ' phenoName ' at location: ' locationPathName]);
end

%
% Store results in a structure to return
%

try 
    regionStruct = struct('DistInRegion',DistInRegion,...
                          'densityInRegion',densityInRegion,...
                          'dataMat',dataMat,...
                          'thetabest',thetabest,...
                          'p0',p0,...
                          'transProbs',transProbs,...
                          'margProbs',margProbs,...
                          'logZ',logZ,...
                          'normEntRegion',normEntRegion,...
                          'ExpectedValSofX',ExpectedValSofX,...
                          'NNcorr',NNcorr);
catch
    disp(['Error in creating structure for phenotype ' phenoName ' at location: ' locationPathName]);
end
