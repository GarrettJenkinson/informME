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
%%%%%%%%%%%             Last Modified: 05/24/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function merges the outputs from MatrixFromBAMfile.m into a single 
% mat file with a single hashtable.
% 
% Example Default Usage:
% MergeMatrices(bamFilename,chr_num)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% MergeMatrices(bamFilename,chr_num,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% bamFilename       Name of the .bam file (without the .bam extension). 
%                   This file must be sorted from least to greatest 
%                   base pair position along the reference seqeunce 
%                   and must be indexed (i.e., the associated .bai file 
%                   must be available).
%
% chr_num           Number representing the chromosome to be processed. 
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MergeMatrices(...,'inputName',inputValue):
%
% pairedEnds     Flag for paired end read support. A value of 1 
%                indicates that the sequencer employed paired end reads, 
%                whereas a value of 0 indicates that the seqeuncer 
%                employed single end reads. Default value: 1.
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               statistical estimation of this phenotype on this
%               chromosome. Default value: 1.
%
% species
%               A string detailing the species from which the data is
%               obtained. e.g. 'Human' or 'Mouse'. Default value: 'Human'.
%
% includeChrInRef
%               A flag telling how the reference chromosomes are named
%                   if 1, then chromosomes are named chr1, chr2, etc. 
%                   if 0, then chromosomes are named 1, 2, etc. 
%               Default value: 0.
%
% CpGlocationPathRoot 
%               Parent path to the files of CpG locations indexed  
%               according to the reference genome in FastaToCpGloc.m.
%               Default Value: './genome/'.
%
% bamFilePathRoot   
%               Parent path to the .bam file. Default value:
%               '../indexedBAMfiles/'.
%
% numBasesToTrim
%               An vector with integer telling how many bases should be 
%               trimmed from the begining of each read. If the vector is of
%               length 2, then the first number tells how many bases to 
%               trim from the first read in a read pair and the second 
%               number tells how many bases should be trimmed from the 
%               second read in the pair. If the vector is of lenght 1 then
%               all reads will have that number of bases trimmed from the 
%               beginning of the read. If no bases should be trimmed, then
%               this should be set equal to 0. Default value: 0.
%
% matricesPathRoot
%               Parent path to write out matrix results. Default value:
%               './matrices/'.
%
% All other optional inputs should only be changed by professionals with a
% detailed understanding of the code.

function MergeMatrices(bamFilename,chr_num,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'bamFilename')
addRequired(p,'chr_num')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'matricesPathRoot',['.' filesep 'matrices' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'includeChrInRef',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','scalar'}))
addParameter(p,'pairedEnds',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','scalar'}))
addParameter(p,'CpGlocationPathRoot',['.' filesep 'genome' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'bamFilePathRoot',['..' filesep 'indexedBAMfiles' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'numBasesToTrim',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'minCpGsReqToModel',10,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))   
       
          
parse(p,bamFilename,chr_num,varargin{:})

species             = p.Results.species;
totalProcessors     = p.Results.totalProcessors;
matricesPathRoot    = p.Results.matricesPathRoot;
includeChrInRef     = p.Results.includeChrInRef;
pairedEnds          = p.Results.pairedEnds;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
bamFilePathRoot     = p.Results.bamFilePathRoot;
numBasesToTrim      = p.Results.numBasesToTrim;
regionSize          = p.Results.regionSize;
minCpGsReqToModel   = p.Results.minCpGsReqToModel;

%
% Manual checks/corrections of inputs
%

if bamFilePathRoot(end)~=filesep
    bamFilePathRoot=[bamFilePathRoot filesep];
end

if matricesPathRoot(end)~=filesep
    matricesPathRoot=[matricesPathRoot filesep];
end

if CpGlocationPathRoot(end)~=filesep
    CpGlocationPathRoot=[CpGlocationPathRoot filesep];
end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all files to verify they exist, create them if not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

chr_num_str = num2str(chr_num);

%
% First check if final result exists. We do this because if accidentally
% invoked twice, this code will serially do all the estimation work that
% was supposed to be done in parallel
%

if exist([matricesPathRoot filesep species filesep 'chr' chr_num_str filesep bamFilename ...
       '.mat'],'file')
    disp('Final merged file already exists.');
    disp('This program will not overwrite an existing file.');
    disp(['In order to recreate this file, first delete existing file: ./matrices' filesep species filesep 'chr' chr_num_str filesep bamFilename ...
       '.mat']);
    return;
end

%
% Check if any of the parallel jobs failed, and if they did, redo the
% computations here
%

for processorNum=1:totalProcessors
    if ~exist([matricesPathRoot species filesep 'chr' chr_num_str filesep bamFilename ...
                num2str(processorNum) '.mat'],'file')
        %
        % File does not exist, redo the computation...must specify all
        % optional inputs just in case user changed one of them from a
        % default value
        %
        MatrixFromBAMfile(bamFilename,chr_num,'CpGlocationPathRoot',CpGlocationPathRoot,'bamFilePathRoot',bamFilePathRoot,...
                          'pairedEnds',pairedEnds,'totalProcessors',totalProcessors,'processorNum',processorNum,...
                           'species',species,'includeChrInRef',includeChrInRef,'numBasesToTrim',numBasesToTrim,...
                           'regionSize',regionSize,'minCpGsReqToModel',minCpGsReqToModel,'matricesPathRoot',matricesPathRoot);
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
% Proceed through all regions in a chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for processorNum = 1:totalProcessors 
    load([matricesPathRoot species filesep 'chr' chr_num_str filesep bamFilename ...
          num2str(processorNum) '.mat'],'mapObjDataTemp');
    mapObjData = [mapObjData;mapObjDataTemp]; %#ok<AGROW>
end

%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save output to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

save([matricesPathRoot  species filesep 'chr' chr_num_str filesep bamFilename ...
       '.mat'],'mapObjData','-v7.3');
  
%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete old files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
  
for processorNum = 1:totalProcessors 
    delete([matricesPathRoot species filesep 'chr' chr_num_str filesep bamFilename ...
          num2str(processorNum) '.mat']);
end  