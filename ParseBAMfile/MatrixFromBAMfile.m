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
% This function processes a .bam file with aligned reads to a reference 
% genome and produces methylation information for nonoverlapping regions 
% of 3000 base pairs long in a given chromosome. The final output for each 
% region is a matrix with -1,0,1 values. Each row of the matrix is a 
% read and each column represents a CpG site within the region. 
% -1 indicates no methylation information is available for the CPG site, 
% 0 indicates that the CpG site is unmethylated, and 1 indicates that the 
% CpG site is methylated. THIS FUNCTION DEPENDS ON A WORKING INSTALLATION
% OF "SAMTOOLS" THAT IS ON THE SYSTEM $PATH.
%
% WARNING: Before running this function, you must run FastaToCpGloc.m ONCE!
%
% Example Default Usage:
% MatrixFromBAMfile(bamFilename,chr_num)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% MatrixFromBAMfile(bamFilename,chr_num,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% bamFilename     - Name of the .bam file (without the .bam extension). 
%                   This file must be sorted from least to greatest 
%                   base pair position along the reference seqeunce 
%                   and must be indexed (i.e., the associated .bai file 
%                   must be available). File name must not contain
%                   unnecessary "." characters, but can contain "_"
%                   instead. File name should end in a character and not a
%                   number, and should be unique from other files. 
%
% chr_num         - Number representing the chromosome to be processed. 
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MergeMatrices(...,'inputName',inputValue):
%
% pairedEnds        Flag for paired end read support. A value of 1 
%                   indicates that the sequencer employed paired end reads, 
%                   whereas a value of 0 indicates that the seqeuncer 
%                   employed single end reads. Default value: 1.
%
% totalProcessors
%               An integer that specifies the number of processors on a
%               computing cluster that will be devoted to the task of
%               statistical estimation of this phenotype on this
%               chromosome. Default value: 1.
%
% processorNum
%               An integer from 1 to totalProcessors which labels which
%               processor of the totalProcessors is being called here. Each
%               processor is assigned a disntinct location of regions
%               within the chromosome to perform estimation on. Default
%               value: 1.
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
% CpGlocationPathRoot 
%               Parent path to the files of CpG locations indexed  
%               according to the reference genome in FastaToCpGloc.m.
%               Default Value: './genome/'.
%
% bamFilePathRoot   
%               Parent path to the .bam file. Default value:
%               '../indexedBAMfiles/'.
%
% matricesPathRoot
%               Parent path to write out matrix results. Default value:
%               './matrices/'.
%
% Any other optional parameters should only be modified by professionals
% with a detailed understanding of the code.

function MatrixFromBAMfile(bamFilename,chr_num,varargin)
                       
                       
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
addParameter(p,'processorNum',1,...
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
processorNum        = p.Results.processorNum;
matricesPathRoot    = p.Results.matricesPathRoot;
includeChrInRef     = p.Results.includeChrInRef;
pairedEnds          = p.Results.pairedEnds;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
bamFilePathRoot     = p.Results.bamFilePathRoot;
numBasesToTrim      = p.Results.numBasesToTrim;
regionSize          = p.Results.regionSize; % base pairs in each region being modeled
minCpGsReqToModel   = p.Results.minCpGsReqToModel; % minimum number of CpG sites in 
                                                   % a region for which methylation 
                                                   % data will be collected
                                                   
%
% Manual checks/corrections of inputs
%

if ~isempty(strfind(bamFilename,'.'))
    disp('ERROR: Input bam file name has period characters.')
    disp('Please rename file without periods. Underscore is allowed.')
    return;
elseif isstrprop(bamFilename(end),'digit')
    disp('ERROR: Input bam file name ends with a number.')
    disp('Please rename file to end with a letter.')
    return;
end

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
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

bamFile           = [bamFilePathRoot species filesep bamFilename '.bam'];
CpGfileNameRoot   = 'CpGlocationChr'; % root file name of .mat files
newLineChar       = sprintf('\n');    % constant used for text processing

%
% get name of chromsome as a string
%
chr_str = num2str(chr_num);


%
% Check if the final output already exists, and exit with error if so
%

if exist([matricesPathRoot species filesep 'chr' chr_str filesep bamFilename ...
           num2str(processorNum) '.mat'],'file')
    disp('Exiting. Final output file already exists:');
    disp([matricesPathRoot species filesep 'chr' chr_str filesep bamFilename ...
      num2str(processorNum) '.mat']);
    return;
elseif exist([matricesPathRoot species filesep 'chr' chr_str filesep  ...
              bamFilename '.mat'],'file')
    disp('Exiting. Final merged output file already exists:');
    disp([matricesPathRoot species filesep 'chr' chr_str filesep  ...
              bamFilename '.mat']);
    return;
end


%
% load chromosome CpG location information
%
load([CpGlocationPathRoot species filesep CpGfileNameRoot chr_str '.mat'],'CpGlocation');

numCpGs = int64(length(CpGlocation));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find which regions are assigned to this processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

finalCpGloc=int64(CpGlocation(numCpGs));

% find all start base-pairs for regions to be modeled on this chromsome
allStartBPs           = int64(1):regionSize:int64(finalCpGloc); % only need to go to the last CpG site, not last BP

% break up the indices of allStartBPs to those relevant to this processor
indexOfProcessor      = processorNum:totalProcessors:length(allStartBPs);

% find all start base-pairs for regions to be modeled on this chromosome by
% this processor
thisProcessorStartBPs = allStartBPs(indexOfProcessor);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize hashtable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

mapObjDataTemp = containers.Map('KeyType','char','ValueType','any');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proceed through all regions in a chromosome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

for startBP = thisProcessorStartBPs 
    try
        %                                only need to go to the last CpG site, 
        %                                not last base pair
        endBP = int64(startBP+regionSize-int64(1)); 
        %                                last region might be too large in 
        %                                base pairs but this does not matter

        %
        % find all CpG sites in the region 
        %
        [lower_CpGindex,upper_CpGindex] = findSortedIndices(CpGlocation,...
                                                            startBP,endBP);
        CpGlocInRegion  = CpGlocation(lower_CpGindex:upper_CpGindex);

        %
        % check if there is a sufficient number of CpG sites in region 
        %
        if length(CpGlocInRegion) < minCpGsReqToModel 
            % skip region since there are not enough CpG sites to model
            continue;
        end

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gather all reads relevant to the current region using samtools
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Define the (possibly smaller than 3000 bp) subregion where reads 
        % will be collected. Reads with overlap to this region may contain 
        % observations of CpG sites within the current 3000 bp genomic region. 
        % We only care about the region having CpG sites (this is why not using
        % startBP & endBP). +1 is assigned to the last CpG site so we do not 
        % miss a read on the reverse complementary strand.
        %
        if includeChrInRef
            region_str = ['chr' chr_str ':' num2str(CpGlocInRegion(1)) '-'...
                          num2str(CpGlocInRegion(end)+1)];
        else
            region_str = [chr_str ':' num2str(CpGlocInRegion(1)) '-' ...
                          num2str(CpGlocInRegion(end)+1)]; 
        end

        %
        % Make the samtools view command with the following options
        % -f 3    makes sure only reads that have both paired ends mapped 
        %         are counted [3 = (2^0)+(2^1)]
        % -q 30   keeps reads with Phred mapping alignment score >= 30 
        %         (i.e., error probability <= 1/1000)
        % -F 3328 excludes PCR duplicates, secondary alignments, and 
        %         chimeric alignments [3328 = (2^8)+(2^10)+(2^11)]
        %
        if pairedEnds 
            command      = ['samtools view -q 30 -F 3328 -f 3 ' bamFile ' ' region_str];
            commandCount = ['samtools view -c -q 30 -F 3328 -f 3 ' bamFile ' ' region_str]; 
        else 
            command      = ['samtools view -q 30 -F 3329 ' bamFile ' ' region_str];
            commandCount = ['samtools view -q 30 -F 3329 ' bamFile ' ' region_str];
            % ignore -f since no paired ends
            % add (2^0) to -F to exclude a paired end read
        end

        %
        % Before running command, check that region is not a overly highly
        % mapped region. If too many reads are mapped to this region, then this
        % indicates that the mapping is unreliable (e.g., due to repeating
        % elements that cannot be uniquely mapped to a reference genome)
        %

        [status,countedNumReadsStr] = system(commandCount);


        if status~=0    % skip this region because something went wrong 
                        % with samtools
            disp(['Error in samtools read at ' region_str]);
            continue;
        else
            countedNumReads=str2double(countedNumReadsStr);
        end
        if countedNumReads<=0
            continue;
        elseif countedNumReads>5000
            disp(['Too many reads mapped to region: ' region_str]);
            continue;
        end

        %
        % Now run actual command, given that the number of reads is normal
        %

        [status,SAMreads] = system(command); % run samtools command

        if status~=0    % skip this region because something went wrong 
                        % with samtools
            disp(['Error in samtools read at ' region_str]);
            continue;
        elseif isempty(SAMreads)
            continue;
        end

        %
        % parse the output of samtools (where each read is separated by a
        % newline character) to have each read as its own cell element.
        % SAMreadsCell will be a 1x1 cell. 
        % SAMreadsCell{1} is column cell with each read in a cell element
        %
        SAMreadsCell = textscan(SAMreads,'%s','delimiter',newLineChar); 

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Process the reads relevant to the current region
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        observedMatrix = MatrixFromReads(SAMreadsCell{1},CpGlocInRegion,pairedEnds,numBasesToTrim); 

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write output hashtable 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        keyStr = ['chr' chr_str '/bp' num2str(startBP) '-' num2str(endBP)];

        dataStruct = struct('observedMatrix',observedMatrix,...
                            'CpGlocInRegion',CpGlocInRegion);

        mapObjDataTemp(keyStr) = dataStruct;                
    
    catch exception
        disp(['Error in MatrixFromBAM in region start bp' num2str(startBP)])
        disp(getReport(exception,'extended','hyperlinks','off'));
    end
end

%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save output to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% check that output directory exists, create if not
if ~exist([matricesPathRoot species filesep 'chr' chr_str],'dir')
    mkdir([matricesPathRoot species filesep 'chr' chr_str]);
end

% write out to file
save([matricesPathRoot species filesep 'chr' chr_str filesep bamFilename ...
      num2str(processorNum) '.mat'],'mapObjDataTemp','-v7.3');
