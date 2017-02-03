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
%%%%%%%%                   MergeEstParams.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 12/01/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function merges the output of EstParamsForChr.m into a single
% hashtable.
%
% USAGE (default):
%
% MergeEstParams(bamFileNames,chr_num,phenoName)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters:
% MergeEstParams(bamFileNames,chr_num,phenoName,'species','Mouse')
%
% MADATORY INPUTS:
%
% bamFileNames
%               A cell array of strings of BAM file names. These BAM files
%               must be processed first by the MatrixFromBAMfile.m function
%               in the ParseBAMfile subdirectory.
%
% chr_num
%               Chromosome number 1 to 22 (in humans) specifying the 
%               chromosome for which statistical estimation must be 
%               performed.
%
% phenoName
%               A string that specifies the name of the phenotype.
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
%               statistical estimation of this phenotype on this
%               chromosome.
%               Default value: 1
%
% estResultsPathRoot
%               A string that specifies the path to the directory in which  
%               the estimation and analysis results are written. 
%               Default value './results/'
%
% genomePathRoot
%               A string that specifies the path to the directory that 
%               contains the results of genome analysis performed by 
%               MatrixFromBamfile.m. 
%               Default value: '../ParseBAMfile/'
%
% regionSize
%               The size of the genomic region used for parameter 
%               estimation (in number of base-pairs).
%               Default value: 3000 
%
% The dafault value of regionSize should only be changed by an expert with 
% a detailed understanding of the code and the methods used. 
%

function MergeEstParams(bamFileNames,chr_num,phenoName,varargin)

% Parse values passed as inputs to the fuction and validate them.

p = inputParser;

addRequired(p,'bamFileNames')
addRequired(p,'chr_num')
addRequired(p,'phenoName')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))
addParameter(p,'estResultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))
          
parse(p,bamFileNames,chr_num,phenoName,varargin{:})

species            = p.Results.species;
totalProcessors    = p.Results.totalProcessors;
estResultsPathRoot = p.Results.estResultsPathRoot;
genomePathRoot     = p.Results.genomePathRoot;
regionSize         = p.Results.regionSize;

% Manual checks/corrections of inputs

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end

if estResultsPathRoot(end)~=filesep
    estResultsPathRoot=[estResultsPathRoot filesep];
end

chr_num_str = num2str(chr_num);

% Loop through all files to verify they exist, create them if not

% First check if final result exists. This is done because, if accidentally
% invoked twice, this function will serially do all estimation work that 
% was supposed to be done in parallel.

if exist([estResultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '.mat'],'file')
    disp('Final merged file already exists.');
    disp('This program will not overwrite an existing file.');
    disp(['In order to recreate this file, first delete existing file: ' ...
          estResultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '.mat']);
    return;
end

% Check if any of the parallel jobs failed, and if they did, redo the
% computations here.

for processorNum = 1:totalProcessors
    if ~exist([estResultsPathRoot species filesep 'chr' ... 
            chr_num_str filesep phenoName '_chr' chr_num_str ...
            '_file' num2str(processorNum) '.mat'],'file')
    
        % File does not exist - redo the computation. Include all optional
        % input values in case user has changed one of the default values.
    
        EstParamsForChr(bamFileNames,chr_num,phenoName,'species',species,...
                       'totalProcessors',totalProcessors,'processorNum',...
                       processorNum,'estResultsPathRoot',...
                       estResultsPathRoot,'genomePathRoot',...
                       genomePathRoot,'regionSize',regionSize);
    end
end

% Initialize hashtable.

mapObjData = containers.Map('KeyType','char','ValueType','any');

% Loop through all files appending results to mapObjData.

for processorNum = 1:totalProcessors
    load([estResultsPathRoot species filesep 'chr' chr_num_str filesep ... 
        phenoName '_chr' chr_num_str '_file' num2str(processorNum) ...
        '.mat'],'mapObjTemp');
    mapObjData = [mapObjData;mapObjTemp]; %#ok<AGROW>
end

% Save output to file.

save([estResultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '.mat'],'mapObjData','-v7.3')
  
% Delete all temporary files.

for processorNum=1:totalProcessors
    delete([estResultsPathRoot species filesep 'chr' chr_num_str filesep ...
        phenoName '_chr' chr_num_str '_file' num2str(processorNum) ...
        '.mat']);
end

end