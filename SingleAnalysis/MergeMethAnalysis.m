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
%%%%%%%%                    MergeMethAnalysis.m                    %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/08/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function merges the output of MethAnalysisForChr.m into a single
% hashtable.
%
% USAGE (default):
%
% MergeMethAnalysis(chr_num,phenoName)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters.
% MergeMethAnalysis(chr_num,phenoName,'species','Mouse')
%
% MANDATORY INPUTS:
%
% chr_num
%               Chromosome number 1 to 22 (in humans) specifying the 
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
%				Default value: 'Human'
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               methylation analysis of this phenotype on this chromosome
%               Default value: 1
%
% resultsPathRoot
%               A string that specifies the path of the directory in which  
%               the methylation analysis results are written. 
%				Default value './results/'
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
%				The ratio regionSize/subregionSize must be an integer.
% 				Default value: 150
%
% The default values of regionSize and subregionSize should only be 
% changed by an expert with a detailed understanding of the code and 
% the methods used. 
%

function MergeMethAnalysis(chr_num,phenoName,varargin)

% Parse values passed as inputs to the fuction and validate them.

p = inputParser;

addRequired(p,'chr_num')
addRequired(p,'phenoName')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},...
					{'nonempty','integer','positive','scalar'}))
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))  
addParameter(p,'estResultsPathRoot',...
               ['..' filesep 'Modeling' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},...
					{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),...
               @(x)validateattributes(x,{'numeric'},...
					{'nonempty','integer','positive','scalar'}))
addParameter(p,'ESIflag',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MCflag',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
           
parse(p,chr_num,phenoName,varargin{:})

species            = p.Results.species;
totalProcessors    = p.Results.totalProcessors;
resultsPathRoot    = p.Results.resultsPathRoot;
estResultsPathRoot = p.Results.estResultsPathRoot;
genomePathRoot     = p.Results.genomePathRoot;
regionSize         = p.Results.regionSize;
subregionSize      = p.Results.subregionSize;
ESIflag            = p.Results.ESIflag;
MCflag             = p.Results.MCflag;

% Manual checks/corrections of inputs.

if resultsPathRoot(end)~=filesep
    resultsPathRoot=[resultsPathRoot filesep];
end

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end

if estResultsPathRoot(end)~=filesep
    estResultsPathRoot=[estResultsPathRoot filesep];
end

chr_num_str = num2str(chr_num);

% Loop through all files to verify they exist, create them if not.

% First check if final result exists. This is done because, if accidentally
% invoked twice, this function will serially do all the estimation work that
% was supposed to be done in parallel.

if exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '_Analysis.mat'],'file')
    disp('Final merged file already exists.');
    disp('This program will not overwrite an existing file.');
    disp(['In order to recreate this file, first delete existing file: ' ...
		  resultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '_Analysis.mat']);
    return;
end

% Check if any of the parallel jobs failed, and if they did, redo the
% computations here.

for processorNum = 1:totalProcessors
    
    if ~exist([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
        phenoName '_Analysis_file' num2str(processorNum) '.mat'],'file')
        
        % File does not exist - redo the computation. Include all optional 
        % values in case user has changed one of the default values.
        
        MethAnalysisForChr(chr_num,phenoName,...
            			    'species',species,...
						  	'totalProcessors',totalProcessors,...
                          	'processorNum',processorNum,...
                            'resultsPathRoot',resultsPathRoot,...
							'estResultsPathRoot',estResultsPathRoot,...
                          	'genomePathRoot',genomePathRoot,...
                            'ESIflag',ESIflag,...
							'MCflag',MCflag,...
							'regionSize',regionSize,...
                          	'subregionSize',subregionSize);
    end
    
end

% Initialize hashtable.

mapObjData = containers.Map('KeyType','char','ValueType','any');

% Loop through all files.

for processorNum = 1:totalProcessors
    load([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
          phenoName '_Analysis_file' num2str(processorNum) '.mat'], ...
          'mapObjTemp');
    mapObjData = [mapObjData;mapObjTemp]; %#ok<AGROW>
end

% Save output to file.

save([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
      phenoName '_Analysis.mat'],'mapObjData','-v7.3')
  
% Delete all temporary files.

for processorNum=1:totalProcessors
    delete([resultsPathRoot species filesep 'chr' chr_num_str filesep ...
            phenoName '_Analysis_file' num2str(processorNum) '.mat']);
end  

end