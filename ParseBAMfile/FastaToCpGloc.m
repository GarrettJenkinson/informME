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
% This function analyzes the reference genome. It finds all CpG sites in the 
% genome and stores their location (we use the convention of a 1-based 
% coordinate system, in which the first base of the sequence is one, and 
% define the location of a CpG site by the position of the C on the forward
% strand of the reference genome). The script also computes the density of 
% CpG sites around each individual site and calculates the distance between
% neighboring CpG sites. This script must be run (only) once before any
% .bam file is processed.
%
% Example Default Usage:
% FastaToCpGloc(FASTAfilename)
%
% Example Usage Modifying optional parameter "species" as name-value pair:
% FastaToCpGloc(FASTAfilename,'species','Mouse')
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% referenceFASTAgenome - Name of FASTA-formatted reference genome to which 
%                        the .bam file has been aligned. 
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MergeMatrices(...,'inputName',inputValue):
%
% species
%               A string detailing the species from which the data is
%               obtained. e.g. 'Human' or 'Mouse'. Default value: 'Human'.
%
% CpGlocationPathRoot 
%               Parent path to the files of CpG locations indexed  
%               according to the reference genome in FastaToCpGloc.m.
%               Default Value: './genome/'.
%
% maxChrNum
%               Maximum number of chromosomes to be processed. The code
%               will process chromosome 1,...,maxChrNum. Defaul Value: 22.
%
% L
%               Window size used in CpG density calculations.
%
% Output to file:
%
% Written to maxChrNum files called CpGlocationChr{chr_num_str}.mat with  
% each containing the vectors CpGlocation, Dist, and density. {chr_num_str} 
% is a string that indicates the chromosome number, 
% whereas CpGlocation is the base-pair location of the C in the CpG 
% site on the reference strand associated with the chromosome {chr_num_str}.
%

function FastaToCpGloc(FASTAfilename,varargin)


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'FASTAfilename')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'CpGlocationPathRoot',['.' filesep 'genome' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'L',1000,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'})) 
addParameter(p,'maxChrNum',22,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'})) 
         
parse(p,FASTAfilename,varargin{:})

species             = p.Results.species;
CpGlocationPathRoot = p.Results.CpGlocationPathRoot;
L                   = p.Results.L;
maxChrNum           = p.Results.maxChrNum;
% species='Human';
%L = 1000;       % window size for density computation
%maxChrNum=22;


%
% Manual checks/corrections of inputs
%

if CpGlocationPathRoot(end)~=filesep
    CpGlocationPathRoot=[CpGlocationPathRoot filesep];
end


%
% Default Human FASTA file downloaded from:
% http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta 
% FASTAfilename = 'Homo_sapiens_assembly19.fasta';
%

referenceFASTAgenome = [CpGlocationPathRoot species filesep FASTAfilename];

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following constants determine how Matlab will represent a string 
% of nucleotides as uint8. Following is not used: A_int = nt2int('A');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
C_int = nt2int('C'); 
G_int = nt2int('G'); 
T_int = nt2int('T');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sequentially proceed through each chromosome in genome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for chr_num=1:maxChrNum
    % convert to string to represent appropriate chromosome
    chr_num_str = num2str(chr_num);
      
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read seqeuence of current chromosome
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    [~,sequenceStr] = fastaread(referenceFASTAgenome,'Blockread',chr_num);
    
    intData = nt2int(sequenceStr); % converts string of all nucleotides 
                                   % on the chromosome to uint8
    clear sequenceStr
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find all CpG sites for the current chromosome
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    prevNC      = T_int;    % this variable indicates the location of the 
                            % previous nucleotide - initialized not to 
                            % be a C
    CpGcount    = int64(0);
    CpGlocation = zeros(1e7,1,'int64'); % initialize locations of CpG sites
    
    chrLength   = int64(length(intData));% compute length of chromosome
    
    for location = int64(1):chrLength
        
        % get location of current nucleotide
        currNC = intData(location);
        
        % check for CpG site
        if prevNC == C_int && currNC== G_int 
            % CpG found
            % previous location was a CpG site (we consider 
            % the C to be the location of the CpG site)
            CpGcount=CpGcount+1;
            CpGlocation(CpGcount) = location-1; 
            % subtract 1 since we are at the G and want location of C 
        end
        prevNC=currNC;    
    end
    CpGlocation = CpGlocation(1:CpGcount);
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make Dist variable for distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    Dist = [CpGlocation(2:end)-CpGlocation(1:(end-1));
            int64(4294967295)]; %#ok<NASGU> 
                                % max value of uint32 ... just chosen as 
                                % a large number since the distance 
                                % to the next CpG site is undefined for  
                                % the last CpG site on a chromosome                              
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute density variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    density = zeros(CpGcount,1);    % preallocate   
    for cpgNum=int64(1):int64(CpGcount)
        density(cpgNum) = nDensity(CpGlocation,cpgNum,L); 
                                    % compute density at current CpG site 
    end
    
    finalCpGloc = CpGlocation(end); %#ok<NASGU>
        
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write output file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    save([CpGlocationPathRoot species filesep 'CpGlocationChr' chr_num_str '.mat'],...
         'CpGlocation','density','Dist','finalCpGloc','chrLength');
    
end