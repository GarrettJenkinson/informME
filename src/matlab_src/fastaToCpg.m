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
%%%%%%%%                     FastaToCpG.m                          %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
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
% fastaToCpg(FASTAfilename)
%
% USAGE (optional): 
%
% Example of optional usage with additional input parameters.  
% FastaToCpG(FASTAfilename,'maxChrNum',23)
%
% MADATORY INPUT:
%
% FASTAfilename
%               Full path of FASTA-formatted reference genome to which  
%               available BAM files have been aligned to.
%
% OPTIONAL INPUTS: 
%
% outdir
%               Path where the output will be stored at. 
% 		Default: "$REFGENEDIR"
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

function fastaToCpg(FASTAfilename,varargin)
% Parse values passed as inputs to the fuction and validate them.
p = inputParser;
addRequired(p,'FASTAfilename')
addParameter(p,'outdir',['.' filesep 'genome' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'maxChrNum',22,@(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'})) 
addParameter(p,'wsize',1000,@(x)validateattributes(x,{'numeric'},{'nonempty',...
               'integer','positive','scalar'})) 

parse(p,FASTAfilename,varargin{:})
outdir              = p.Results.outdir;
maxChrNum           = p.Results.maxChrNum;
wsize               = p.Results.wsize;

% The following constants determine how MATLAB will represent a string 
% of nucleotides as uint8 (only C, G and T nucleotides are used).
C_int = nt2int('C'); 
G_int = nt2int('G'); 
T_int = nt2int('T');

% Sequentially proceed through each chromosome in the genome.
for chr_num = 1:maxChrNum
    chr_num_str = num2str(chr_num); % Convert to string to represent 
                                    % appropriate chromosome.
                                    
    [~,sequenceStr] = fastaread(FASTAfilename,'Blockread',chr_num); 
                                    % Read seqeuence of current chromosome.

    intData = nt2int(sequenceStr); % Convert string of nucleotides 
                                   % in the chromosome to uint8.
    clear sequenceStr

    % Find all CpG sites for the current chromosome.
    % prevNC variable indicates the location of the 
    % previous nucleotide - initialized not to  be a C.
    prevNC      = T_int;    

    CpGcount    = int64(0);
    CpGlocation = zeros(1e7,1,'int64');   % Initialize locations of CpG sites.
    chrLength   = int64(length(intData)); % Compute length of chromosome.
    for location = int64(1):chrLength
	% Get location of current nucleotide.
        currNC = intData(location);         
        % Check for CpG site.
        if prevNC == C_int && currNC == G_int 
            % CpG found. Previous location was a CpG site  
            % (consider C o be location of the CpG site).
            CpGcount = CpGcount+1;
            % Subtract 1 since at G and want location of C. 
            CpGlocation(CpGcount) = location-1; 
        end
        prevNC = currNC; 
    end
    CpGlocation = CpGlocation(1:CpGcount);
    
    % Compute CpG distances.
    Dist = [CpGlocation(2:end)-CpGlocation(1:(end-1));
            int64(4294967295)]; 
	%#ok<NASGU> 
        % Max value of uint32. Just chosen as 
        % a large number since the distance 
        % to the next CpG site is undefined for  
        % the last CpG site in a chromosome.                               
   
    % Compute CpG densities.
    density = zeros(CpGcount,1);
    for cpgNum=int64(1):int64(CpGcount)
       	% Compute density at current CpG site. 
        density(cpgNum) = nDensity(CpGlocation,cpgNum,wsize); 
    end
    finalCpGloc = CpGlocation(end); %#ok<NASGU>
        
    % Save results. 
    save([outdir filesep 'CpGlocationChr' chr_num_str '.mat'],'CpGlocation','density', ...
	  'Dist','finalCpGloc','chrLength');
    
end

