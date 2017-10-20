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
%%%%%%%%                   mergeEstParams.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 12/01/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function merges the output of EstParamsForChr.m into a single
% hashtable.
%
% USAGE (default):
%
% mergeEstParams(matrices_path,reference_path,estimation_path,chr_num,prefix)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters:
% mergeEstParams(matrices_path,reference_path,estimation_path,chr_num,prefix,...
%	'outdir','/path/to/output')
%
% MADATORY INPUTS:
%
% matrices_path
%               A string that specifies the path to the directory in which  
%               the matrices are stored. 
%               Default value "$INTERMEDIATE"
%
% reference_path
%               A string that specifies the path to the directory that 
%               contains the results of genome analysis performed by 
%               MatrixFromBamfile.m. 
%               Default value: "$REFGENEDIR"
%
% estimation_path
%               A string that specifies the path to the directory in which  
%               the estimation and analysis results are written. 
%               Default value "$INTERMEDIATE"
%
% chr_num
%               Chromosome number 1 to 22 (in humans) specifying the 
%               chromosome for which statistical estimation must be 
%               performed.
%
% prefix
%               A string that specifies the name of the phenotype.
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
% regionSize
%               The size of the genomic region used for parameter 
%               estimation (in number of base-pairs).
%               Default value: 3000 
%
% The dafault value of regionSize should only be changed by an expert with 
% a detailed understanding of the code and the methods used. 
%

function mergeEstParams(matrices_path,reference_path,estimation_path,chr_num,prefix,varargin)

% Parse values passed as inputs to the fuction and validate them.
p = inputParser;
addRequired(p,'matrices_path')
addRequired(p,'reference_path')
addRequired(p,'estimation_path')
addRequired(p,'chr_num')
addRequired(p,'prefix')
addParameter(p,'outdir',['.' filesep 'results' filesep],@(x)validateattributes(x,...
		{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,@(x)validateattributes(x,{'numeric'},{'nonempty',...
               	'integer','positive','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},{'nonempty',...
               	'integer','positive','scalar'}))
          
parse(p,matrices_path,reference_path,estimation_path,chr_num,prefix,varargin{:})
totalProcessors     = p.Results.totalProcessors;
regionSize          = p.Results.regionSize;
outdir              = p.Results.outdir;

% Manual checks/corrections of inputs
if matrices_path(end)~=filesep
    matrices_path=[matrices_path filesep];
end
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if estimation_path(end)~=filesep
    estimation_path=[estimation_path filesep];
end
if outdir(end)~=filesep
    outdir=[outdir filesep];
end

chr_num_str = num2str(chr_num);

% Loop through all files to verify they exist, create them if not
% First check if final result exists. This is done because, if accidentally
% invoked twice, this function will serially do all estimation work that 
% was supposed to be done in parallel.
outfile=[outdir 'chr' chr_num_str filesep prefix '_fit.mat'];
if exist(outfile,'file')
    fprintf(2,'Final merged file already exists.\n');
    fprintf(2,'This program will not overwrite an existing file.\n');
    fprintf(2,['In order to recreate this file, first delete existing file: ' ...
          outdir 'chr' chr_num_str filesep prefix '_fit.mat\n']);
    return;
end

% Check if any of the parallel jobs failed, and if they did, redo the
% computations here.
for processorNum = 1:totalProcessors
    temp_file=[estimation_path 'chr' chr_num_str filesep ...
	    prefix '_params_chr' chr_num_str '_file' num2str(processorNum) '.mat'];
    if ~exist(temp_file,'file')
        % File does not exist - redo the computation. Include all optional
        % input values in case user has changed one of the default values.    
        fprintf(2,'WARNING: The following file does not exist and is being processed:');
        fprintf(2,[temp_file '\n']);
        estParamsForChr(matrices_path,reference_path,chr_num,prefix,'totalProcessors',totalProcessors,...
                        'processorNum',processorNum,'outdir',estimation_path,'regionSize',regionSize);
    end
end

% Initialize hashtable.
mapObjData = containers.Map('KeyType','char','ValueType','any');

% Loop through all files appending results to mapObjData.
for processorNum = 1:totalProcessors
    load([estimation_path 'chr' chr_num_str filesep prefix '_params_chr' ...
        chr_num_str '_file' num2str(processorNum) '.mat'],'mapObjTemp');
    mapObjData = [mapObjData;mapObjTemp]; %#ok<AGROW>
end

% Save output to file.
save([outdir 'chr' chr_num_str filesep prefix '_fit.mat'],'mapObjData','-v7.3')

% Delete all temporary files.
for processorNum=1:totalProcessors
    delete([estimation_path 'chr' chr_num_str filesep prefix '_params_chr' ...
            chr_num_str '_file' num2str(processorNum) '.mat']);
end

