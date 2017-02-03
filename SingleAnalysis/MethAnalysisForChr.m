% informME: An information-theoretic pipeline for WGBS data
% Copyright (C) 2017, Garrett Jenkinson (jenkinson@jhu.edu)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                   MethAnalysisForChr.m                    %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/08/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function performs methylation analysis of a given chromosome in a 
% single phenotype. The function can be used on a computing cluster to 
% break the analysis work to many independent parallel job processes. 
% This is performed only after EstParamsForChr.m in the Modeling 
% subdirectory is run to build the Ising models for the phenotype.
%
% USAGE (default):
%
% MethAnalysisForChr(chr_num,phenoName)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters.
% MethAnalysisForChr(chr_num,phenoName,'species','Mouse')
%
% MANDATORY INPUTS:
%
% chr_num
%               Chromosome number (1 to 22 in humans) specifying the 
%               chromosome for which methylation analyis must be 
%               performed.
%
% phenoName
%               A string that specifies the name of the phenotype to be 
%               analyzed. 
%
% OPTIONAL INPUTS: 
%
% species
%               A string that specifies the species from which the data is
%               obtained (e.g., 'Human' or 'Mouse'). 
%               Default value: 'Human'
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               methylation analysis of this phenotype on this chromosome
%               Default value: 1
%
% processorNum
%               An integer from 1 to totalProcessors which labels which
%               processor of the totalProcessors is being called. Each
%               processor is assigned a distinct location of regions
%               within the chromosome to perform methylation analysis. 
%               Default value: 1
%
% resultsPathRoot
%               A string that specifies the path of the directory in which  
%               the methylation analysis results are written. 
%               Default value './results/'
%
% estResultsPathRoot
%               A string that specifies the path to the directory that 
%               contains the results of parameter estimation performed 
%               by EstParamsForChr.m.
%               Default value '../Modeling/results/'
%
% genomePathRoot
%               A string that specifies the path to the directory that 
%               contains the results of analysis of the reference genome 
%               performed by FastaToCpG.m as well as the results of 
%               methylation calling performed by MatrixFromBAMfile.m.
%               Default value: '../ParseBAMfile/'
%
% ESIflag
%               Flag that determines whether this function performs 
%               computation of the entropic sensitivity index (ESI). 
%               0: no ESI computation. 
%               1: allow ESI computation.
%               Default value: 1
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
% regionSize
%               The size of the genomic regions used for parameter 
%               estimation (in number of base pairs).
%               Default value: 3000 
%
% subregionSize
%               The size of the subregions of a genomic region used 
%               for methylation analysis (in number of base pairs). 
%               The ratio regionSize/subregionSize must be an integer.
%               Default value: 150
%
% The default values of regionSize and subregionSize should only be 
% changed by an expert with a detailed understanding of the code and 
% the methods used. 
%

function MethAnalysisForChr(chr_num,phenoName,varargin)

% Parse values passed as inputs to the fuction and validate them.

p = inputParser;

addRequired(p,'chr_num')
addRequired(p,'phenoName')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},...
				{'nonempty','integer','positive','scalar'}))
addParameter(p,'processorNum',1,...
               @(x)validateattributes(x,{'numeric'},...
				{'nonempty','integer','positive','scalar'}))
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'estResultsPathRoot',...
				['..' filesep 'Modeling' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'ESIflag',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MCflag',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'regionSize',int64(3000),... 
               @(x)validateattributes(x,{'numeric'},...
				{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),... 
               @(x)validateattributes(x,{'numeric'},...
				{'nonempty','integer','positive','scalar'}))

parse(p,chr_num,phenoName,varargin{:})

species            = p.Results.species;
totalProcessors    = p.Results.totalProcessors;
processorNum       = p.Results.processorNum;
resultsPathRoot    = p.Results.resultsPathRoot;
estResultsPathRoot = p.Results.estResultsPathRoot;
genomePathRoot     = p.Results.genomePathRoot;
ESIflag            = p.Results.ESIflag;
MCflag             = p.Results.MCflag;
regionSize         = p.Results.regionSize;
subregionSize      = p.Results.subregionSize;

% Manual checks/corrections of inputs.

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end

if resultsPathRoot(end)~=filesep
    resultsPathRoot=[resultsPathRoot filesep];
end

if estResultsPathRoot(end)~=filesep
    estResultsPathRoot=[estResultsPathRoot filesep];
end

if processorNum > totalProcessors
    disp('error: processorNum must be <= totalProcessors');
    return;
end

% Convert chr_num to string.

chr_num_str=num2str(chr_num);

% Check if the final output already exists, and exit with error if so.

if exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis_file' num2str(processorNum) '.mat'],'file')
    disp('Exiting. Final output file already exists:');
    disp([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis_file' num2str(processorNum) '.mat']);
    return;
elseif exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis.mat'],'file')
    disp('Exiting. Final merged output file already exists:');
    disp([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis.mat']);
    return;
end

% Find last CpG site on chromosome. 

CpGdata = [genomePathRoot 'genome' filesep species filesep ...
			'CpGlocationChr' chr_num_str '.mat'];

load(CpGdata,'finalCpGloc','Dist','density','CpGlocation');

% Find which regions are assigned to this processor.

% Find all start base pairs for regions to be modeled on this chromosome.
% Only need to go to the last CpG site, not last base pair.
allStartBPs = int64(1):regionSize:int64(finalCpGloc); 

% Break up the indices of allStartBPs to those relevant to this processor.
indexOfProcessor = processorNum:totalProcessors:length(allStartBPs);

% Find all start base pairs for regions to be modeled on this chromosome 
% by this processor.
thisProcessorStartBPs = allStartBPs(indexOfProcessor);

% Load data for this chromosome

load([estResultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '.mat'],'mapObjData');

% Loop through all regions on chromosome designated to this processor.

% Initialize hashtable
mapObjTemp = containers.Map('KeyType','char','ValueType','any');

for startBP = thisProcessorStartBPs
    try
		% The last region might be too large in base pairs but this 
        % does not matter.
        endBP = int64(startBP+regionSize-int64(1));

        % Find relative path of this region.
        locationPathName = ['chr' chr_num_str '/bp' ...
                               num2str(startBP) '-' num2str(endBP)];
        
        if isKey(mapObjData,locationPathName)
            % Find CpG sites. 
            [lower_index,upper_index] = ...
                       findSortedIndices(CpGlocation,startBP,endBP);
            CpGlocs_local = CpGlocation(lower_index:upper_index);

            localEstStruct = mapObjData(locationPathName);

            % Estimate parameters for this region.
            FDregionStruct = MethAnalysisForRegion(localEstStruct,CpGlocs_local,...
							   startBP,endBP,subregionSize,ESIflag,MCflag);
            mapObjTemp(locationPathName) = FDregionStruct;
        end
    catch ME
        disp(['Error in MethAnalysisForRegion at ' locationPathName ' with stack:']) 
        ME.stack.file
        ME.stack.name
        ME.stack.line
        ME.identifier
    end
end

% Check that results folder exists, create if not.

if ~exist([resultsPathRoot species filesep 'chr' chr_num_str],'dir')
    mkdir([resultsPathRoot species filesep 'chr' chr_num_str]);
end

% Save output to file.

save([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis_file' num2str(processorNum) '.mat'],...
      'mapObjTemp','-v7.3')
  
end
