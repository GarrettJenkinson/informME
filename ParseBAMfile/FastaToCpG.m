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
%%%%%%%%                     FastaToCpG.m                          %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 11/30/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is used to analyze a reference genome in order to find 
% and store the locations of all CpG sites within each chromosome and 
% to compute the CpG densities at each CpG site as well as the distances 
% between neighboring CpG sites. A 1-based coordinate system is used, 
% in which the first base is assigned to position 1 and the location 
% of a CpG site is defined by the position of the C nucleotide on the 
% forward strand of the reference genome. 
%
% This function must be run ONLY ONCE before proceeding with analysis  
% of BAM files.
%
% USAGE (default):
%
% FastaToCpG(FASTAfilename)
%
% USAGE (optional): 
%
% Example of optional usage with additional input parameters.  
% FastaToCpG(FASTAfilename,'species','Mouse')
%
% MADATORY INPUT:
%
% FASTAfilename
%               Name of a FASTA-formatted reference genome to which  
%               available BAM files have been aligned. This file must 
%               be placed in the subdirectory determined by the optional 
%               input CpGlocationPathRoot (see below).
%
% OPTIONAL INPUTS: 
%
% species
%               A string that specifies the species from which the data is 
%               obtained (e.g., 'Human' or 'Mouse'). 
%               Default value: 'Human'
%
% CpGlocationPathRoot 
%               Path to the root subdirectory where the FASTA files for 
%               the reference genomes are located and where the outputs 
%               of this function are stored.
%               Default value: './genome/'
%
% maxChrNum
%               Maximum number of chromosomes to be processed. The  
%               function will process chromosomes 1,2,...,maxChrNum. 
%               Default value: 22
%
% wsize
%               Window size used in CpG density calculations.
%               Default value: 1000
%

function FastaToCpG(FASTAfilename,varargin)

% Parse values passed as inputs to the fuction and validate them.

p = inputParser;

addRequired(p,'FASTAfilename')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'CpGlocationPathRoot',['.' filesep 'genome' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'maxChrNum',22,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'})) 
addParameter(p,'wsize',1000,...
               @(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'})) 

parse(p,FASTAfilename,varargin{:})

species             = p.Results.species;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
maxChrNum           = p.Results.maxChrNum;
wsize               = p.Results.wsize;

% Manual check/correction of inputs.

if CpGlocationPathRoot(end)~=filesep
    CpGlocationPathRoot=[CpGlocationPathRoot filesep];
end

% Load FASTA file.

referenceFASTAgenome = [CpGlocationPathRoot species filesep FASTAfilename];

% The following constants determine how MATLAB will represent a string 
% of nucleotides as uint8 (only C, G and T nucleotides are used).

C_int = nt2int('C'); 
G_int = nt2int('G'); 
T_int = nt2int('T');

% Sequentially proceed through each chromosome in the genome.

for chr_num = 1:maxChrNum
    
    chr_num_str = num2str(chr_num); % Convert to string to represent 
                                    % appropriate chromosome.
                                    
    [~,sequenceStr] = fastaread(referenceFASTAgenome,'Blockread',chr_num); 
                                    % Read seqeuence of current chromosome.

    intData = nt2int(sequenceStr); % Convert string of nucleotides 
                                   % in the chromosome to uint8.
    clear sequenceStr
    
    % Find all CpG sites for the current chromosome.
    
    prevNC      = T_int;    % This variable indicates the location of the 
                            % previous nucleotide - initialized not to 
                            % be a C.
    CpGcount    = int64(0);
    CpGlocation = zeros(1e7,1,'int64');   % Initialize locations of CpG sites.
    chrLength   = int64(length(intData)); % Compute length of chromosome.
    
    for location = int64(1):chrLength
        
        currNC = intData(location); % Get location of current nucleotide.
        
        % Check for CpG site.
        
        if prevNC == C_int && currNC == G_int 
                            % CpG found. Previous location was a CpG site  
                            % (consider C o be location of the CpG site).
            CpGcount = CpGcount+1;
            CpGlocation(CpGcount) = location-1; 
                            % Subtract 1 since at G and want location of C. 
        end
        
        prevNC = currNC; 
        
    end
    
    CpGlocation = CpGlocation(1:CpGcount);
    
    % Compute CpG distances.
    
    Dist = [CpGlocation(2:end)-CpGlocation(1:(end-1));
            int64(4294967295)]; %#ok<NASGU> 
                                % Max value of uint32. Just chosen as 
                                % a large number since the distance 
                                % to the next CpG site is undefined for  
                                % the last CpG site in a chromosome.                               
   
    % Compute CpG densities.
    
    density = zeros(CpGcount,1);
    
    for cpgNum=int64(1):int64(CpGcount)
        density(cpgNum) = nDensity(CpGlocation,cpgNum,wsize); 
                                   % Compute density at current CpG site. 
    end
    
    finalCpGloc = CpGlocation(end); %#ok<NASGU>
        
    % Save results. 
   
    save([CpGlocationPathRoot species filesep 'CpGlocationChr'...
        chr_num_str '.mat'],'CpGlocation','density','Dist',...
        'finalCpGloc','chrLength');
    
end
