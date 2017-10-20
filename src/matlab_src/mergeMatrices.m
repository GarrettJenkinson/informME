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
%%%%%%%%                    mergeMatrices.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/30/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function merges the outputs from MatrixFromBAMfile.m into a single 
% MATLAB MAT file comprised of a single hashtable.
% 
% USAGE (default):
% mergeMatrices(bam_prefix,chr_num)
%
% USAGE (optional):
% Example of optional usage with additional input parameters.
% mergeMatrices(bam_prefix,chr_num,'CpGlocationPathRoot','/path/to/ref')
%
% MADATORY INPUTS:
%
% bam_prefix       
%                Name of the BAM file (without the .bam extension). This 
%                file must be sorted from the least to the greatest base 
%                pair position along the reference sequence and must be 
%                indexed (i.e., the associated BAI file must be available). 
%                The file name must not contain "." characters, but can 
%                contain "_" instead. Moreover, the file name should be 
%                unique.
%
% chr_num        
%                Number representing the chromosome to be processed. 
%
% OPTIONAL INPUTS:
%
% totalProcessors
%                An integer that specifies the number of processors on a
%                computing cluster that will be devoted to the task of
%                building data matrices for this choromosome. 
%                Default value: 1
%
% CpGlocationPathRoot 
%                Parent path to the files of CpG locations indexed  
%                according to the reference genome in FastaToCpGloc.m.
%                Default value: "$REFGENEDIR"
%
% bamFilePathRoot   
%                Parent path to the subdirectory where the BAM file is 
%		 located.
%                Default value: "$BAMDIR"
%
% matricesPathRoot
%                Parent path to the subdirectory where the output of this  
%                function is stored. 
%                Default value: "$INTERMEDIATE"
%
% pairedEnds     
%                Flag for paired end read support. A value of 1 indicates 
%                that the sequencer employed paired end reads, whereas a 
%                value of 0 indicates that the sequencer employed single 
%                end reads. 
%                Default value: 1
%
% numBasesToTrim
%                A vector of integers specifying how many bases should be 
%                trimmed from the begining of each read. If the vector 
%                contains two integers, then the first integer specifies 
%                how many bases to trim from the first read in a read pair, 
%                whereas the second integer specifies how many bases should 
%                be trimmed from the second read in the pair. If the 
%                vector contains one integer, then all reads will have 
%                that number of bases trimmed from the beginning of the 
%                read. If no bases are to be trimmed, then this input 
%                must be set to 0. 
%                Default value: 0
%
% includeChrInRef
%                A flag specifying how the reference chromosomes are named
%                   if 1, then chromosomes are named chr1, chr2, etc. 
%                   if 0, then chromosomes are named 1, 2, etc. 
%                Default value: 0
%
% regionSize     
%                The size of the genomic regions for which methylation 
%                information is produced (in number of base pairs).
%                Default value: 3000
%
% minCpGsReqToModel
%                The minimum number of CpG sites within a genomic region 
%                required to perform statistical estimation.
%                Default value: 10
%
% The default values of regionSize and minCpGsReqToModel should only be 
% changed by an expert with a detailed understanding of the code and the 
% methods used. 

function mergeMatrices(bam_prefix,chr_num,varargin)

% Parse values passed as inputs to the fuction and validate them.
p = inputParser;
addRequired(p,'bam_prefix')
addRequired(p,'chr_num')
addParameter(p,'totalProcessors',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'CpGlocationPathRoot',['.' filesep 'genome' filesep],...
                @(x)validateattributes(x,{'char'},{'nonempty'}))   
addParameter(p,'outdir',['.' filesep 'genome' filesep],...
                @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'bamFilePathRoot',['..' filesep 'indexedBAMfiles' filesep],...
                @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'matricesPathRoot',['.' filesep 'matrices' filesep],...
                @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'pairedEnds',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','scalar'}))       
addParameter(p,'numBasesToTrim',0,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer'}))  
addParameter(p,'includeChrInRef',0,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,...
		{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'minCpGsReqToModel',10,@(x)validateattributes(x,...
		{'numeric'},{'nonempty','integer','positive','scalar'}))
                              
parse(p,bam_prefix,chr_num,varargin{:})

totalProcessors     = p.Results.totalProcessors;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
bamFilePathRoot     = p.Results.bamFilePathRoot;
matricesPathRoot    = p.Results.matricesPathRoot;
outdir              = p.Results.outdir;
pairedEnds          = p.Results.pairedEnds;
numBasesToTrim      = p.Results.numBasesToTrim;
includeChrInRef     = p.Results.includeChrInRef;
regionSize          = p.Results.regionSize;
minCpGsReqToModel   = p.Results.minCpGsReqToModel;

% Manual checks/corrections of inputs.
if bamFilePathRoot(end)~=filesep
    bamFilePathRoot=[bamFilePathRoot filesep];
end
if matricesPathRoot(end)~=filesep
    matricesPathRoot=[matricesPathRoot filesep];
end
if CpGlocationPathRoot(end)~=filesep
    CpGlocationPathRoot=[CpGlocationPathRoot filesep];
end

% Get name of chromsome as a string.
chr_num_str = num2str(chr_num);

% Loop through all files to verify they exist, create them if not.
% First, check if final result exists. If it does, do not recompute.
if exist([outdir 'chr' chr_num_str filesep bam_prefix '_matrices.mat'],'file')
    fprintf(2,'Final merged file already exists.\n');
    fprintf(2,'This program will not overwrite the existing file.\n');
    fprintf(2,['To recreate this file, first delete existing file:' ...
        outdir 'chr' chr_num_str filesep bam_prefix '_matrices.mat\n']);
    return;
end

% Check if any of the parallel jobs failed, and if they did, redo the
% computations here.
for processorNum=1:totalProcessors
    processor_file=[matricesPathRoot 'chr' chr_num_str ...
            filesep bam_prefix '_matrices_file_' num2str(processorNum) '.mat'];
    if ~exist(processor_file,'file')
        % File does not exist, redo the computation. Must specify all
        % optional inputs just in case user changed one of them from a
        % default value.
        fprintf(2,['WARNING: The following file did not exist and is being ' ... 
                'created: ' processor_file '\n']);
        matrixFromBam(bam_prefix,chr_num,...
		      'CpGlocationPathRoot',CpGlocationPathRoot,...
		      'bamFilePathRoot',bamFilePathRoot,...
                      'pairedEnds',pairedEnds,...
                      'totalProcessors',totalProcessors,...
                      'processorNum',processorNum,...
                      'outdir',matricesPathRoot,...
                      'includeChrInRef',includeChrInRef,...
                      'numBasesToTrim',numBasesToTrim,...
                      'regionSize',regionSize,...
                      'minCpGsReqToModel',minCpGsReqToModel)
    end
end

% Initialize hashtable.
mapObjData = containers.Map('KeyType','char','ValueType','any');

% Proceed through all regions in a chromosome. 
for processorNum = 1:totalProcessors 
    load([matricesPathRoot 'chr' chr_num_str filesep ...
        bam_prefix '_matrices_file_' num2str(processorNum) '.mat'],'mapObjDataTemp');
    mapObjData = [mapObjData;mapObjDataTemp]; %#ok<AGROW>
end

% Save output to file. 
save([outdir filesep 'chr' chr_num_str filesep ...
        bam_prefix '_matrices.mat'],'mapObjData','-v7.3');
  
% Delete old files. 
for processorNum = 1:totalProcessors 
    delete([matricesPathRoot 'chr' chr_num_str filesep ...
        bam_prefix '_matrices_file_' num2str(processorNum) '.mat']);
end  

