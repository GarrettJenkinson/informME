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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  informME: Information-Theoretic Analysis of Methylation  %%%%%%%%
%%%%%%%%                   MatrixFromBAMfile.m                     %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/30/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function processes a BAM file with aligned reads to a reference 
% genome and produces methylation information for nonoverlapping genomic 
% regions (containing the same number of base pairs) in a given chromosome. 
% The final output for each genomic region is a matrix with -1,0,1 values. 
% Each row of the matrix is a methylation read, whereas each column 
% represents a CpG site within the genomic region. A value of -1 indicates 
% no methylation information is available for the CPG site, 0 indicates 
% that the CpG site is unmethylated, and 1 indicates that the CpG site 
% is methylated. 
% 
% This function depends on a working instalation of SAMtools that is on  
% the system path $PATH.
%
% Before running this function, FastaToCpG.m must be run ONCE. 
%
% USAGE (default):
%
% MatrixFromBAMfile(bamFilename,chr_num)
%
% USAGE (optional): 
%
% Example of optional usage with additional input parameters.
% MatrixFromBAMfile(bamFilename,chr_num,'species','Mouse')
%
% MADATORY INPUTS:
%
% bamFilename
%                Name of the BAM file (without the .bam extension). This 
%                file must be sorted from the least to the greatest base 
%                pair position along the reference sequence and must be 
%                indexed (i.e., the associated BAI file must be available). 
%                The file name must not contain "." characters, but can 
%                contain "_" instead. Moreover, the file name should be 
%                unique and end with a letter (not a number). 
%
% chr_num
%                Number representing the chromosome to be processed. 
%
% OPTIONAL INPUTS:
%
% species
%                A string that specifies the species from which the data is 
%                obtained (e.g., 'Human' or 'Mouse'). 
%                Default value: 'Human'
%
% totalProcessors
%                An integer that specifies the number of processors on a
%                computing cluster that will be devoted to the task of
%                building data matrices for this choromosome. 
%                Default value: 1
%
% processorNum
%                An integer from 1 to totalProcessors that labels which
%                processor is being called. Each processor is assigned 
%                distinct locations of genomic regions within the 
%                chromosome that are used to generate data matrices. 
%                Default value: 1
%
% CpGlocationPathRoot 
%                Path to the root subdirectory where the outputs of this 
%                function are stored.
%                Default value: './genome/'
%
% bamFilePathRoot   
%                Path to the subdirectory where the BAM file is located.
%                Default value: './indexedBAMfiles/'
%
% matricesPathRoot
%                Path to the subdirectory where the output of this function 
%                is stored. 
%                Default value: './matrices/'
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

function MatrixFromBAMfile(bamFilename,chr_num,varargin)
                                           
% Parse values passed as inputs to the fuction and validate them.

p = inputParser;

addRequired(p,'bamFilename')
addRequired(p,'chr_num')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'totalProcessors',1,...
               @(x)validateattributes(x,{'numeric'},...
               {'nonempty','integer','positive','scalar'}))
addParameter(p,'processorNum',1,...
               @(x)validateattributes(x,{'numeric'},...
               {'nonempty','integer','positive','scalar'}))    
addParameter(p,'CpGlocationPathRoot',['.' filesep 'genome' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))     
addParameter(p,'bamFilePathRoot',['..' filesep 'indexedBAMfiles' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))                
addParameter(p,'matricesPathRoot',['.' filesep 'matrices' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'pairedEnds',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','scalar'}))         
addParameter(p,'numBasesToTrim',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer'}))              
addParameter(p,'includeChrInRef',0,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','scalar'}))              
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))     
addParameter(p,'minCpGsReqToModel',10,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'}))           
          
parse(p,bamFilename,chr_num,varargin{:})

species             = p.Results.species;
totalProcessors     = p.Results.totalProcessors;
processorNum        = p.Results.processorNum;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
bamFilePathRoot     = p.Results.bamFilePathRoot;
matricesPathRoot    = p.Results.matricesPathRoot;
pairedEnds          = p.Results.pairedEnds;
numBasesToTrim      = p.Results.numBasesToTrim;
includeChrInRef     = p.Results.includeChrInRef;
regionSize          = p.Results.regionSize;       
minCpGsReqToModel   = p.Results.minCpGsReqToModel;
                                                   
% Manual checks/corrections of inputs. 

if ~isempty(strfind(bamFilename,'.'))
    disp('ERROR: Input BAM file name contains period characters.')
    disp('Please rename file without periods. Underscore is allowed.')
    return;
elseif isstrprop(bamFilename(end),'digit')
    disp('ERROR: Input BAM file name ends with a number.')
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
           
% Initialize. 

bamFile           = [bamFilePathRoot species filesep bamFilename '.bam'];
CpGfileNameRoot   = 'CpGlocationChr';
newLineChar       = sprintf('\n');     % Constant used for text processing. 

% Get name of chromsome as a string.

chr_num_str = num2str(chr_num);

% Check if final output already exists, and exit with error if so.

if exist([matricesPathRoot species filesep 'chr' chr_num_str filesep ...
        bamFilename num2str(processorNum) '.mat'],'file')
    disp('Exiting. Final output file already exists:');
    disp([matricesPathRoot species filesep 'chr' chr_num_str filesep ...
        bamFilename num2str(processorNum) '.mat']);
    return;
elseif exist([matricesPathRoot species filesep 'chr' chr_num_str filesep ...
              bamFilename '.mat'],'file')
    disp('Exiting. Final merged output file already exists:');
    disp([matricesPathRoot species filesep 'chr' chr_num_str filesep ...
              bamFilename '.mat']);
    return;
end

% Load chromosome CpG location information.

load([CpGlocationPathRoot species filesep CpGfileNameRoot chr_num_str ...
    '.mat'],'CpGlocation');

numCpGs = int64(length(CpGlocation));

% Find which regions are assigned to this processor. 

finalCpGloc=int64(CpGlocation(numCpGs));

% Find all start base pairs for regions to be modeled on this chromsome.

allStartBPs = int64(1):regionSize:int64(finalCpGloc); 
                % Only need to go to the last CpG site, not last base pair.

% Break up the indices of allStartBPs to those relevant to this processor.

indexOfProcessor = processorNum:totalProcessors:length(allStartBPs);

% Find all start base pairs for regions to be modeled on this chromosome by
% this processor.

thisProcessorStartBPs = allStartBPs(indexOfProcessor);

% Initialize hashtable.

mapObjDataTemp = containers.Map('KeyType','char','ValueType','any');

% Proceed through all regions in a chromosome. 

for startBP = thisProcessorStartBPs 
    try
                                    % Only need to go to the last CpG site, 
                                    % not last base pair. 
        endBP = int64(startBP+regionSize-int64(1)); 
                                    % Last region might be too large in 
                                    % base pairs but this does not matter. 
                                    
        [lower_CpGindex,upper_CpGindex] = findSortedIndices(CpGlocation,...
                                                            startBP,endBP);
                                    % Find all CpG sites in the region.                                          
        CpGlocInRegion = CpGlocation(lower_CpGindex:upper_CpGindex);

        % Check if there is a sufficient number of CpG sites in region. 
        
        if length(CpGlocInRegion) < minCpGsReqToModel 
            % Skip region since there are not enough CpG sites to model. 
            continue;
        end

        % Gather all reads relevant to the current genomic region using 
        % SAMtools. Determine the subregion where reads will be collected. 
        % Reads with overlap to this subregion may contain observations of 
        % CpG sites within the current genomic region. We only care about 
        % the subregion having CpG sites (this is why not using startBP 
        % & endBP). +1 is assigned to the last CpG site so we do not 
        % miss a read on the reverse complementary strand.
        
        if includeChrInRef
            region_str = ['chr' chr_num_str ':' num2str(CpGlocInRegion(1)) '-'...
                          num2str(CpGlocInRegion(end)+1)];
        else
            region_str = [chr_num_str ':' num2str(CpGlocInRegion(1)) '-' ...
                          num2str(CpGlocInRegion(end)+1)]; 
        end

        % Make the SAMtools view command with the following options: 
        % -f 3    makes sure only reads that have both paired ends mapped 
        %         are counted [3 = (2^0)+(2^1)].
        % -q 30   keeps reads with Phred mapping alignment score >= 30 
        %         (i.e., error probability <= 1/1000).
        % -F 3328 excludes PCR duplicates, secondary alignments, and 
        %         chimeric alignments [3328 = (2^8)+(2^10)+(2^11)].
        
        if pairedEnds 
            command      = ['samtools view -q 30 -F 3328 -f 3 ' ...
                             bamFile ' ' region_str];
            commandCount = ['samtools view -c -q 30 -F 3328 -f 3 ' ...
                             bamFile ' ' region_str]; 
        else 
            command      = ['samtools view -q 30 -F 3329 ' ...
                             bamFile ' ' region_str];
            commandCount = ['samtools view -q 30 -F 3329 ' ...
                             bamFile ' ' region_str];
                            % Ignore -f since no paired ends.
                            % Add (2^0) to -F to exclude a paired end read.
        end

        % Before running SAMtools view command, check whether too many 
        % reads are mapped to the genomic region. This will indicate 
        % that the mapping is unreliable (e.g., due to repetitive 
        % elements that cannot be uniquely mapped to a reference genome).

        [status,countedNumReadsStr] = system(commandCount);
        
        if status~=0      % Skip this region because something went wrong  
                          % with SAMtools.
            disp(['Error in SAMtools read at ' region_str]);
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

        % Now run actual SAMtools view command, given that the number of 
        % reads is normal.

        [status,SAMreads] = system(command);

        if status~=0       % Skip this region because something went wrong 
                           % with SAMtools.
            disp(['Error in SAMtools read at ' region_str]);
            continue;
        elseif isempty(SAMreads)
            continue;
        end

        % Parse the output of SAMtools (where each read is separated by a
        % newline character) to have each read as its own cell element.
        % SAMreadsCell will be a 1x1 cell. 
        % SAMreadsCell{1} is column cell with each read in a cell element.
      
        SAMreadsCell = textscan(SAMreads,'%s','delimiter',newLineChar); 

        
        % Process the reads relevant to the current region.
        observedMatrix = MatrixFromReads(SAMreadsCell{1},...
                              CpGlocInRegion,pairedEnds,numBasesToTrim); 

        % Write output hashtable. 
        keyStr = ['chr' chr_num_str '/bp' num2str(startBP) '-' num2str(endBP)];
        dataStruct = struct('observedMatrix',observedMatrix,...
                            'CpGlocInRegion',CpGlocInRegion);
        mapObjDataTemp(keyStr) = dataStruct;                
    
    catch exception
        disp(['Error in MatrixFromBAM in region start bp' num2str(startBP)])
        disp(getReport(exception,'extended','hyperlinks','off'));
    end
end

% Check that output directory exists, create if not.

if ~exist([matricesPathRoot species filesep 'chr' chr_num_str],'dir')
    mkdir([matricesPathRoot species filesep 'chr' chr_num_str]);
end

% Write output to file.

save([matricesPathRoot species filesep 'chr' chr_num_str filesep bamFilename ...
      num2str(processorNum) '.mat'],'mapObjDataTemp','-v7.3');
