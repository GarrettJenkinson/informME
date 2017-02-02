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
% EstParamsForChr(bamFileNames,chr_num,phenoName)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% EstParamsForChr(bamFileNames,chr_num,phenoName,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% bamFileNames
%               A cell array of strings of BAM file names. These BAM files
%               should have already been processed by MatrixFromBAMfile.m
%               to produce matrices. e.g. {'bamFilenameWithoutExtension'}
%
% chr_num
%               The chromosome number 1 to 22 (in humans) which
%               details the chromosome for which statistical estimation
%               should be performed.
%
% phenoName
%               A string which details the name of the phenotype that is
%               represented by the BAM files specified in bamFileNames
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% EstParamsForChr(...,'inputName',inputValue):
%
% species
%               A string detailing the species from which the data is
%               obtained. i.e., 'Human' or 'Mouse'. Default: 'Human'
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               statistical estimation of this phenotype on this chromosome
%               Default value:1
%
% processorNum
%               An integer from 1 to totalProcessors which labels which
%               processor of the totalProcessors is being called here. Each
%               processor is assigned a disntinct location of regions
%               within the chromosome to perform estimation on. Default
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

function EstParamsForChr(bamFileNames,chr_num,phenoName,varargin)%species,totalProcessors,processorNum)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'bamFileNames')
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

     
parse(p,bamFileNames,chr_num,phenoName,varargin{:})

species         = p.Results.species;
totalProcessors = p.Results.totalProcessors;
processorNum    = p.Results.processorNum;
resultsPathRoot = p.Results.resultsPathRoot;
genomePathRoot  = p.Results.genomePathRoot;
regionSize      = p.Results.regionSize;


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




chr_num_str = num2str(chr_num);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the final output already exists, and exit with error if so
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_chr' chr_num_str '_file' num2str(processorNum) '.mat'],'file')
    disp('Exiting. Final output file already exists:');
    disp([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_chr' chr_num_str '_file' num2str(processorNum) '.mat']);
    return;
elseif exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
                phenoName '_chr' chr_num_str '.mat'],'file')
    disp('Exiting. Final merged output file already exists:');
    disp([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
                phenoName '_chr' chr_num_str '.mat']);
    return;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find which regions are assigned to this processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Load CpG site data on chromosome (to determine how many regions to loop through, etc.)
CpGdata = [genomePathRoot 'genome' filesep species filesep 'CpGlocationChr' chr_num_str '.mat'];
load(CpGdata,'finalCpGloc','Dist','density','CpGlocation');

% find all start base-pairs for regions to be modeled on this chromsome
allStartBPs  = int64(1):regionSize:int64(finalCpGloc); % only need to go to the last CpG site, not last BP

% break up the indices of allStartBPs to those relevant to this processor
indexOfProcessor = processorNum:totalProcessors:length(allStartBPs);

% find all start base-pairs for regions to be modeled on this chromosome by
% this processor
thisProcessorStartBPs = allStartBPs(indexOfProcessor);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data for this chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for dataFile = 1:length(bamFileNames)
    load([genomePathRoot 'matrices' filesep species filesep 'chr' chr_num_str ... 
          filesep bamFileNames{dataFile} '.mat'],'mapObjData');
    dataMapObj{dataFile} = mapObjData; %#ok<AGROW>
    clear mapObjData;
end

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
        
        %
        % load data for region
        %
        dataMat=[];
        for dataFile = 1:length(bamFileNames)
            if isKey(dataMapObj{dataFile},locationPathName)
                regionDataStruct = dataMapObj{dataFile}(locationPathName);
                if isempty(dataMat) %get values of Dist and density in region 
                    [lower_index,upper_index] = findSortedIndices(CpGlocation,startBP,endBP);
                    DistInRegion    = Dist(lower_index:upper_index);
                    densityInRegion = density(lower_index:upper_index);
                end
                dataMat = [dataMat;regionDataStruct.observedMatrix]; %#ok<AGROW>
            end
        end
        
        %
        % Estimate parameters
        %
        if ~isempty(dataMat)
            % estimate parameters for this region
            regionStruct = EstimateParams(locationPathName,phenoName,DistInRegion,densityInRegion,dataMat);

            % Add to hashtable
            if ~isempty(regionStruct)
                mapObjTemp(locationPathName) = regionStruct;
            end
        end
    catch
        disp(['Error in EstimateParams at ' locationPathName]) 
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save output to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% check that results folder exists, create it if not
if ~exist([resultsPathRoot species filesep 'chr' chr_num_str],'dir')
    mkdir([resultsPathRoot species filesep 'chr' chr_num_str]);
end

% write output to file
save([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_chr' chr_num_str '_file' num2str(processorNum) '.mat'],'mapObjTemp','-v7.3')
  
end