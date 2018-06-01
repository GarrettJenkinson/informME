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
%%%%%%%%  informME: Information-Theoretic analysis of Methylation  %%%%%%%%
%%%%%%%%                   methAnalysisForChr.m                    %%%%%%%%
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
% methAnalysisForChr(prefix,chr_num,reference_path,estimation_path)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters.
% methAnalysisForChr(prefix,chr_num,reference_path,estimation_path,
% 		     'outdir','/path/to/output')
%
% MANDATORY INPUTS:
%
% prefix
%               A string that specifies the name of the phenotype to be 
%               analyzed. 
%
% chr_num
%               Chromosome number (1 to 22 in humans) specifying the 
%               chromosome for which methylation analyis must be 
%               performed.
%
% reference_path
%               A string that specifies the path to the directory that 
%               contains the results of analysis of the reference genome 
%               performed by FastaToCpG.m as well as the results of 
%               methylation calling performed by MatrixFromBAMfile.m.
%		Default: "$REFGENEDIR"
%
% estimation_path
%               A string that specifies the path to the directory that 
%               contains the results of parameter estimation performed 
%               by EstParamsForChr.m.
%		Default: "$INTERMEDIATE"
%
% OPTIONAL INPUTS: 
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
% outdir
%               A string that specifies the path of the directory in which  
%               the methylation analysis results are written. 
%               Default value './results/'
%
% ESIflag
%               Flag that determines whether this function performs 
%               computation of the entropic sensitivity index (ESI). 
%               0: no ESI computation. 
%               1: allow ESI computation.
%               Default value: 1
%
% MSIflag
%               Flag that determines whether this function performs
%               computation of the methylation sensitivity index (MSI).
%               0: no MSI computation.
%               1: allow MSI computation.
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

function methanalysisForChr(prefix,chr_num,reference_path,estimation_path,varargin)

% Parse values passed as inputs to the fuction and validate them.
p = inputParser;
addRequired(p,'prefix')
addRequired(p,'chr_num')
addRequired(p,'reference_path')
addRequired(p,'estimation_path')
addParameter(p,'outdir',['.' filesep 'results' filesep],...
               	@(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'ESIflag',0,...
               	@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MSIflag',0,...
               	@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MCflag',0,...
               	@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'totalProcessors',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'processorNum',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))

parse(p,prefix,chr_num,reference_path,estimation_path,varargin{:})

totalProcessors     = p.Results.totalProcessors;
processorNum        = p.Results.processorNum;
outdir              = p.Results.outdir;
ESIflag             = p.Results.ESIflag;
MSIflag             = p.Results.MSIflag;
MCflag              = p.Results.MCflag;
regionSize          = p.Results.regionSize;
subregionSize       = p.Results.subregionSize;

% Manual checks/corrections of inputs.
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if estimation_path(end)~=filesep
    estimation_path=[estimation_path filesep];
end
if outdir(end)~=filesep
    outdir=[outdir filesep];
end
if processorNum > totalProcessors
    fprintf(2,'error: processorNum must be <= totalProcessors\n');
    return;
end

% Convert chr_num to string.
chr_num_str=num2str(chr_num);

% Check if the final output already exists, and exit with error if so.
temp_file=[outdir 'chr' chr_num_str filesep prefix '_analysis_file' ...
            num2str(processorNum) '.mat'];
merged_file=[outdir 'chr' chr_num_str filesep prefix '_analysis.mat'];
if exist(temp_file,'file')
    fprintf(2,'Exiting. Final output file already exists:\n%s\n',temp_file);
    return;
elseif exist(merged_file,'file')
    fprintf(2,'Exiting. Final merged output file already exists:\n%s\n',merged_file);
    return;
end

% Find last CpG site on chromosome. 
CpGdata = [reference_path filesep 'CpGlocationChr' chr_num_str '.mat'];
load(CpGdata,'finalCpGloc','Dist','density','CpGlocation');

% Find all start base pairs for regions to be modeled on this chromosome.
% Only need to go to the last CpG site, not last base pair.
allStartBPs = int64(1):regionSize:int64(finalCpGloc); 

% Break up the indices of allStartBPs to those relevant to this processor.
indexOfProcessor = processorNum:totalProcessors:length(allStartBPs);

% Find all start base pairs for regions to be modeled on this chromosome 
% by this processor.
thisProcessorStartBPs = allStartBPs(indexOfProcessor);

% Load data for this chromosome
load([estimation_path 'chr' chr_num_str filesep prefix '_fit.mat'],'mapObjData');

% Loop through all regions on chromosome designated to this processor.
% Initialize hashtable
mapObjTemp = containers.Map('KeyType','char','ValueType','any');

for startBP = thisProcessorStartBPs
    try
	% The last region might be too large in base pairs but this 
        % does not matter.
        endBP = int64(startBP+regionSize-int64(1));

        % Find relative path of this region.
        locationPathName = ['chr' chr_num_str '/bp' num2str(startBP) '-' num2str(endBP)];
        
        if isKey(mapObjData,locationPathName)
            % Find CpG sites. 
            [lower_index,upper_index] = findSortedIndices(CpGlocation,startBP,endBP);
            CpGlocs_local = CpGlocation(lower_index:upper_index);
            localEstStruct = mapObjData(locationPathName);

            % Estimate parameters for this region.
            FDregionStruct = methAnalysisForRegion(localEstStruct,CpGlocs_local,...
						startBP,endBP,subregionSize,ESIflag,MSIflag,MCflag);
            mapObjTemp(locationPathName) = FDregionStruct;
        end
    catch ME
        fprintf(2,['Error in methanalysisForRegion.m at ' locationPathName ' with message:' ME.message '\n']) 
    end
end

% Check that results folder exists, create if not.
if ~exist([outdir filesep 'chr' chr_num_str],'dir')
    mkdir([outdir filesep 'chr' chr_num_str]);
end

% Save output to file.
save([outdir 'chr' chr_num_str filesep prefix '_analysis_file' ...
    num2str(processorNum) '.mat'],'mapObjTemp','-v7.3')
  
