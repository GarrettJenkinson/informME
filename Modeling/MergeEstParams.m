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
% This function merges the output of EstParamsForChr.m into a single
% hashtable.
%
% Example Default Usage:
% MergeEstParams(chr_num,phenoName,bamFileNames)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% MergeEstParams(chr_num,phenoName,bamFileNames,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
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
% bamFileNames
%               A cell array of strings of BAM file names. These BAM files
%               should have already been processed by MatrixFromBAMfile.m
%               to produce matrices. e.g. {'bamFilenameWithoutExtension'}
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MergeEstParams(...,'inputName',inputValue):
%
% species
%               A string detailing the species from which the data is
%               obtained. e.g. 'Human' or 'Mouse'
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               statistical estimation of this phenotype on this chromosome
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

function MergeEstParams(chr_num,phenoName,bamFileNames,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'chr_num')
addRequired(p,'phenoName')
addRequired(p,'bamFileNames')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
          
          
parse(p,chr_num,phenoName,bamFileNames,varargin{:})

species         = p.Results.species;
totalProcessors = p.Results.totalProcessors;
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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all files to verify they exist, create them if not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

chr_num_str=num2str(chr_num);

%
% First check if final result exists. We do this because if accidentally
% invoked twice, this code will serially do all the estimation work that
% was supposed to be done in parallel
%

if exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '.mat'],'file')
    disp('Final merged file already exists.');
    disp('This program will not overwrite an existing file.');
    disp(['In order to recreate this file, first delete existing file: ' ...
          resultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '.mat']);
    return;
end

%
% Check if any of the parallel jobs failed, and if they did, redo the
% computations here
%

for processorNum=1:totalProcessors
    if ~exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
               phenoName '_chr' chr_num_str '_file' num2str(processorNum) '.mat'],'file')
        %
        % File does not exist, redo the computation. Include all optional
        % input values in case user has changed one of the default values.
        %
        EstParamsForChr(bamFileNames,chr_num,phenoName,'species',species,...
                       'totalProcessors',totalProcessors,'processorNum',processorNum,...
                       'resultsPathRoot',resultsPathRoot,'genomePathRoot',genomePathRoot,...
                       'regionSize',regionSize);
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize hashtable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mapObjData = containers.Map('KeyType','char','ValueType','any');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all files appending results to mapObjData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for processorNum=1:totalProcessors
    load([resultsPathRoot species filesep 'chr' chr_num_str filesep phenoName ...
          '_chr' chr_num_str '_file' num2str(processorNum) '.mat'],'mapObjTemp');
    mapObjData = [mapObjData;mapObjTemp]; %#ok<AGROW>
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save output to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

save([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '.mat'],'mapObjData','-v7.3')
  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete all temporary files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for processorNum=1:totalProcessors
    delete([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
            phenoName '_chr' chr_num_str '_file' num2str(processorNum) '.mat']);
end  

end