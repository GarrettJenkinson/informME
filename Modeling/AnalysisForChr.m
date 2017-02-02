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
%%%%%%%%%%%             Last Modified: 05/22/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function takes a list of bam files (all corresponding to the same
% phenotype) and performs statistical model estimation for regions within
% the specified chromosome of interest. In particular, this function can be
% used on a computing cluster to break up the work to many independent
% parallel job processes. This is performed only after MatrixFromBamfile.m
% has been run to produce the data formats required for statistical
% esimtation.
%
% Example Default Usage:
% AnalysisForChr(chr_num,phenoName)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% AnalysisForChr(chr_num,phenoName,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% chr_num
%               A number with the chromosome for which analysis
%               should be performed.
%
% phenoName
%               A string which details the name of the phenotype on which
%               estimation has already been performed using
%               EstParamsForChr.m
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% AnalysisForChr(...,'inputName',inputValue):
%
% species
%               A string detailing the species from which the data is
%               obtained. i.e., 'Human' or 'Mouse'. Default: 'Human'
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               analysis of this phenotype on this chromosome
%               Default value:1
%
% processorNum
%               An integer from 1 to totalProcessors which labels which
%               processor of the totalProcessors is being called here. Each
%               processor is assigned a disntinct location of regions
%               within the chromosome to perform analysis on. Default
%               value:1
%
% resultsPathRoot
%               A string with the path to the results directory where the
%               estimation and analysis results are written in its
%               subdirectories. Default value './results/'
%
% genomePathRoot
%               A string with the path to the results directory where the
%               genome analysis results are written in its subdirectories.
%               Default value: '../ParseBAMfile/'
%
% Any other optional inputs should only be modified by professionals with a
% detailed understanding of the code.

function AnalysisForChr(chr_num,phenoName,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'chr_num')
addRequired(p,'phenoName')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'processorNum',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'subRegionSize',int64(150),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))

parse(p,chr_num,phenoName,varargin{:})

species         = p.Results.species;
totalProcessors = p.Results.totalProcessors;
processorNum    = p.Results.processorNum;
resultsPathRoot = p.Results.resultsPathRoot;
genomePathRoot  = p.Results.genomePathRoot;
regionSize      = p.Results.regionSize;
subRegionSize   = p.Results.subRegionSize;


%
% Manual checks/corrections of inputs
%

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end
if resultsPathRoot(end)~=filesep
    resultsPathRoot=[resultsPathRoot filesep];
end
if processorNum > totalProcessors
    disp('error: processorNum must be <= totalProcessors');
    return;
end

% convert chr_num to string
chr_num_str=num2str(chr_num);



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the final output already exists, and exit with error if so
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find last CpG site on chromosome (to determine how many regions to loop through)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

CpGdata = [genomePathRoot 'genome' filesep species filesep 'CpGlocationChr' chr_num_str '.mat'];

load(CpGdata,'finalCpGloc','Dist','density','CpGlocation');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find which regions are assigned to this processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% find all start base-pairs for regions to be modeled on this chromsome
allStartBPs           = int64(1):regionSize:int64(finalCpGloc); % only need to go to the last CpG site, not last BP

% break up the indices of allStartBPs to those relevant to this processor
indexOfProcessor      = processorNum:totalProcessors:length(allStartBPs);

% find all start base-pairs for regions to be modeled on this chromosome by
% this processor
thisProcessorStartBPs = allStartBPs(indexOfProcessor);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data for this chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

load([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '.mat'],'mapObjData');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all regions on chromosome designated to this processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%initialize hashtable
mapObjTemp = containers.Map('KeyType','char','ValueType','any');

for startBP = thisProcessorStartBPs
    try
        endBP = int64(startBP+regionSize-int64(1)); % the last region might be too large in BP but this does not matter

        % find relative path of this region
        locationPathName = ['chr' chr_num_str '/bp' num2str(startBP) '-' num2str(endBP)];
        
        if isKey(mapObjData,locationPathName)
            % find CpG sites 
            [lower_index,upper_index] = findSortedIndices(CpGlocation,startBP,endBP);
            CpGlocs_local = CpGlocation(lower_index:upper_index);

            localEstStruct = mapObjData(locationPathName);

            % estimate parameters for this region
            FDregionStruct = AnalysisRegion(localEstStruct,CpGlocs_local,startBP,endBP,subRegionSize);
            mapObjTemp(locationPathName) = FDregionStruct;
        end
    catch ME
        display(['Error in AnalysisRegion at ' locationPathName ' with stack:']) 
        ME.stack.file
        ME.stack.name
        ME.stack.line
        ME.identifier
    end
end

%
% Save output to file
%

save([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis_file' num2str(processorNum) '.mat'],'mapObjTemp','-v7.3')
  
end
