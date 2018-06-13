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
%%%%%%%%                    mergeMethAnalysis.m                    %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function merges the output of MethAnalysisForChr.m into a single
% hashtable.
%
% USAGE (default):
%
% mergeMethAnalysis(analysis_path,prefix,chr_num,reference_path,estimation_path)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters.
% mergeMethAnalysis(analysis_path,prefix,chr_num,reference_path,estimation_path,
%	'outdir','/path/to/output')
%
% MANDATORY INPUTS:
%
% analysis_path
%               A string that specifies the path of the directory in which  
%               the methylation analysis results are written. 
%		Default: "$INTERMEDIATE"
%
% prefix
%               A string that specifies the name of the phenotype to be 
%               analyzed. 
%
% chr_num
%               Chromosome number 1 to 22 (in humans) specifying the 
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
% ESIflag
%               Flag that determines whether this function performs 
%               computation of the entropic sensitivity index (ESI). 
%               0: no ESI computation. 
%               1: allow ESI computation.
%               Default value: 0
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
%		The ratio regionSize/subregionSize must be an integer.
% 		Default value: 150
%
% The default values of regionSize and subregionSize should only be 
% changed by an expert with a detailed understanding of the code and 
% the methods used. 
%

function mergeMethAnalysis(analysis_path,prefix,chr_num,reference_path,estimation_path,varargin)

% Argument parser
p = inputParser;

% Mandatory arguments
addRequired(p,'analysis_path')
addRequired(p,'prefix')
addRequired(p,'chr_num')
addRequired(p,'reference_path')
addRequired(p,'estimation_path')

% Optional arguments
addParameter(p,'outdir',['.' filesep 'merged_analysis' filesep],...
               	@(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'ESIflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty',...
		'scalar'}))
addParameter(p,'MSIflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty',...
		'scalar'}))
addParameter(p,'MCflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty',...
		'scalar'}))
           
% Parse optional arguments
parse(p,analysis_path,prefix,chr_num,reference_path,estimation_path,varargin{:})
outdir              = p.Results.outdir;
totalProcessors     = p.Results.totalProcessors;
regionSize          = p.Results.regionSize;
subregionSize       = p.Results.subregionSize;
ESIflag             = p.Results.ESIflag;
MSIflag             = p.Results.MSIflag;
MCflag              = p.Results.MCflag;

% Manual checks/corrections of inputs.
if analysis_path(end)~=filesep
    analysis_path=[analysis_path filesep];
end
if outdir(end)~=filesep
    outdir=[outdir filesep];
end
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if estimation_path(end)~=filesep
    estimation_path=[estimation_path filesep];
end

chr_num_str = num2str(chr_num);

% Loop through all files to verify they exist, create them if not.
% First check if final result exists. This is done because, if accidentally
% invoked twice, this function will serially do all the estimation work that
% was supposed to be done in parallel.
merged_file=[outdir 'chr' chr_num_str filesep prefix '_analysis.mat'];
if exist(merged_file,'file')
    fprintf(2,'Final merged file already exists.');
    fprintf(2,'This program will not overwrite an existing file.\n');
    fprintf(2,['In order to recreate this file, first delete existing file: ' merged_file '\n']);
    return;
end

% Check if any of the parallel jobs failed, and if they did, redo the
% computations here.
for processorNum = 1:totalProcessors
    temp_file=[analysis_path 'chr' chr_num_str filesep prefix '_analysis_file' ...
                num2str(processorNum) '.mat'];
    if ~exist(temp_file,'file')
        % Spit out which file was missing 
        fprintf(2,['WARNING: The following file is being re-processed: ' temp_file '\n']);

        % File does not exist - redo the computation. Include all optional 
        % values in case user has changed one of the default values.
        methAnalysisForChr(prefix,chr_num,reference_path,estimation_path,...
                           'outdir',analysis_path,'totalProcessors',totalProcessors,...
                           'processorNum',processorNum,'ESIflag',ESIflag,'MSIflag',MSIflag,...
			   'MCflag',MCflag,'regionSize',regionSize,'subregionSize',subregionSize);
    end
    
end

% Initialize hashtable.
mapObjData = containers.Map('KeyType','char','ValueType','any');

% Loop through all files.
for processorNum = 1:totalProcessors
    load([analysis_path 'chr' chr_num_str filesep prefix '_analysis_file' ...
        num2str(processorNum) '.mat'],'mapObjTemp');
    mapObjData = [mapObjData;mapObjTemp]; %#ok<AGROW>
end

% Save output to file.

save([outdir 'chr' chr_num_str filesep prefix '_analysis.mat'],'mapObjData','-v7.3')
  
% Delete all temporary files.
for processorNum=1:totalProcessors
    delete([analysis_path 'chr' chr_num_str filesep prefix '_analysis_file' ...
            num2str(processorNum) '.mat']);
end  

