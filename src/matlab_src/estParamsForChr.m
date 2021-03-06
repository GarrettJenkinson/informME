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
%%%%%%%%                   EstParamsForChr.m                       %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a list of BAM files (which correspond to the same
% phenotype) and performs statistical model estimation within a specific  
% chromosome of interest. The function can be used on a computing cluster 
% to break the work of model estimation to many independent parallel job 
% processes. This is performed only after matrixFromBam.m.
%
% USAGE (default):
%
% estParamsForChr(mat_files,prefix,matrices_path,reference_path,chr_num)
%
% USAGE (optional):
% 
% Example of optional usage with additional input parameters.
% estParamsForChr(mat_files,prefix,matrices_path,reference_path,chr_num,...
%		  'regionSize',2000)
%
% MANDATORY INPUTS:
%
% mat_files
%		All the .mat files to be included in the model. This can be a
% 		single .mat file or multiple files in the form of a comma-sepa-
%		rated list of files.
%
% prefix
%               A string that specifies the name of the modeled phenotype.
%		The output files produced will contain this prefix.
%
% matrices_path
%               A string that specifies the path to the directory that 
%               where the output will be stored.
%		Default: "$INTERMEDIATE"
%
% reference_path
%               A string that specifies the path to the directory that 
%               contains the results of analysis of the reference genome 
%               performed by FastaToCpG.m as well as the results of 
%               methylation calling performed by MatrixFromBAMfile.m.
%		Default: "$REFGENEDIR"
%
% chr_num
%               Chromosome number 1 to 22 (in humans) specifying the 
%               chromosome for which statistical estimation must be 
%               performed.
%		Default: 22.
%
%
% OPTIONAL INPUTS: 
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               statistical estimation of this phenotype on this
%               chromosome.
%               Default value: 1
%
% processorNum
%               An integer from 1 to totalProcessors that labels which
%               processor of the totalProcessors is being called. Each
%               processor is assigned a distinct location of regions
%               within the chromosome to perform estimation. 
%               Default value: 1
%
% regionSize
%               The size of the genomic region used for parameter 
%               estimation (in number of base pairs).
%               Default value: 3000 
%
% The default value of regionSize should only be changed by an expert with 
% a detailed understanding of the code and the methods used. 
%

function estParamsForChr(mat_files,prefix,matrices_path,reference_path,chr_num,varargin)
% Parse values passed as inputs to the fuction and validate them.
p = inputParser;
addRequired(p,'mat_files')
addRequired(p,'prefix')
addRequired(p,'matrices_path')
addRequired(p,'reference_path')
addRequired(p,'chr_num')
addParameter(p,'totalProcessors',1,@(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))
addParameter(p,'processorNum',1,@(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))
addParameter(p,'outdir',['.' filesep 'results' filesep],@(x)validateattributes(x,...
		{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))

parse(p,mat_files,prefix,matrices_path,reference_path,chr_num,varargin{:})

totalProcessors    = p.Results.totalProcessors;
processorNum       = p.Results.processorNum;
estResultsPathRoot = p.Results.outdir;
regionSize         = p.Results.regionSize;

% Manual checks/corrections of inputs.
if matrices_path(end)~=filesep
    matrices_path=[matrices_path filesep];
end
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if estResultsPathRoot(end)~=filesep
    estResultsPathRoot=[estResultsPathRoot filesep];
end
if processorNum > totalProcessors
    fprintf(2,'error: processorNum must be <= totalProcessors\n');
    return;
end

chr_num_str = num2str(chr_num);

% Check if final output already exists, and exit with error if so.
if exist([estResultsPathRoot filesep 'chr' chr_num_str filesep ...
      prefix '_params_chr' chr_num_str '_file' num2str(processorNum) ...
      '.mat'],'file')
    fprintf(2,'Exiting. Final output file already exists:');
    fprintf(2,[estResultsPathRoot filesep 'chr' chr_num_str filesep ...
      prefix '_params_chr' chr_num_str '_file' num2str(processorNum) '.mat\n']);
    return;
elseif exist([estResultsPathRoot filesep 'chr' chr_num_str  ...
               filesep prefix '_fit.mat'],'file')
    fprintf(2,'Exiting. Final merged output file already exists:');
    fprintf(2,[estResultsPathRoot filesep 'chr' chr_num_str filesep ...
                prefix '_fit.mat\n']);
    return;
end

% Find which genomic regions are assigned to this processor.
% Load CpG site data on chromosome (to determine how many regions 
% to loop through, etc.).
CpGdata = [reference_path filesep 'CpGlocationChr' chr_num_str '.mat'];
load(CpGdata,'finalCpGloc','Dist','density','CpGlocation');

% Find all start base-pairs for regions to be modeled on this chromsome.
allStartBPs  = int64(1):regionSize:int64(finalCpGloc); 
               % Only need to go to the last CpG site, not last base pair.

% Break up the indices of allStartBPs to those relevant to this processor.
indexOfProcessor = processorNum:totalProcessors:length(allStartBPs);

% Find all start base pairs for regions to be modeled on this chromosome
% by this processor.
thisProcessorStartBPs = allStartBPs(indexOfProcessor);

% Load data for this chromosome for all matrices provided.
mat_files = strsplit(mat_files,',');
for index = 1:length(mat_files)
    base = char(mat_files(index));
    file_path = [matrices_path 'chr' chr_num_str filesep base '_matrices.mat']
    load(file_path,'mapObjData');
    dataMapObj{index} = mapObjData;
    clear mapObjData;
end

% Loop through all regions on chromosome designated to this processor.
% Initialize hashtable.
mapObjTemp = containers.Map('KeyType','char','ValueType','any');

for startBP = thisProcessorStartBPs
    try
        endBP = int64(startBP+regionSize-int64(1)); 
        % The last region might be too large in base pairs but this 
        % does not matter.

        % Find relative path of this region.
        locationPathName = ['chr' chr_num_str '/bp' ...
            num2str(startBP) '-' num2str(endBP)];
        
        % Load data for region.
        dataMat=[];
	for i = 1:length(mat_files)
            if isKey(dataMapObj{i},locationPathName)
                regionDataStruct = dataMapObj{i}(locationPathName);
                if isempty(dataMat) % Get values of Dist and density in region.  
                    [lower_index,upper_index] = ...
                        findSortedIndices(CpGlocation,startBP,endBP);
                    DistInRegion    = Dist(lower_index:upper_index);
                    densityInRegion = density(lower_index:upper_index);
                end
                % recreate observed matrix from sparse storage structure
                observedMatrix = full(regionDataStruct.sparseObsMatPlusOne)-1;
                dataMat = [dataMat;observedMatrix];%#ok<AGROW>
            end
        end
        
        % Estimate parameters.
        if ~isempty(dataMat)
            % Estimate parameters for this region.
            regionStruct = estimateParams(locationPathName,prefix,...
                DistInRegion,densityInRegion,dataMat);
            % Add to hashtable.
            if ~isempty(regionStruct)
                mapObjTemp(locationPathName) = regionStruct;
            end
        end
    catch ME
        fprintf(2,['Error in estimateParams at ' locationPathName '\n']) 
        fprintf(2,getReport(ME));
    end
end

% Save output to file.
% Check if results folder exists, create it if not.
if ~exist([estResultsPathRoot filesep 'chr' chr_num_str],'dir')
    mkdir([estResultsPathRoot filesep 'chr' chr_num_str]);
end

% Write output to file.
save([estResultsPathRoot filesep 'chr' chr_num_str filesep prefix '_params_chr' ...
      chr_num_str '_file' num2str(processorNum) '.mat'],'mapObjTemp','-v7.3')
  
