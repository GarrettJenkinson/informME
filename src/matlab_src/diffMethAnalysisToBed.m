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
%%%%%%%%               makeBEDsForDiffMethAnalysis.m               %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 06/01/2018                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function makes BED files for the differential version of the 
% methylation analysis results obtained by means of MethAnalysisForChr.m 
% applied on two dinstict phenotypes. 
%
% USAGE (default):
%
% makeBedsForDiffMethAnalysis(prefix_1,prefix_2,analysis_path_1,...
%	analysis_path_2,reference_path)
%
% MANDATORY INPUTS:
%
% prefix_X
%               Strings with the first string specifying 
%               the name of the first phenotype and the second string 
%               specifying the name of the second phenotype used for  
%               differential methylation analysis. Both phenotypes
%               must have already been analyzed with methAnalysisForChr.m.
% 
% analysis_path_X
%               A string that specifies the path of the directory in which  
%               the methylation analysis results obtained by 
%               MethAnalysisForChr.m are stored. 
%		Default: "$INTERMEDIATE"
%
% reference_path
%               A string that specifies the path to the directory that 
%               contains the results of analysis of the reference genome 
%               performed by FastaToCpG.m as well as the results of 
%               methylation calling performed by MatrixFromBAMfile.m.
%		Default: "$REFGENEDIR"
%
% OPTIONAL INPUTS: 
%
% minChrNum
%               A number specifying the starting chromosome that will be 
%               included in the BED files.
%               Default value: 1
%
% maxChrNum
%               A number specifying the last chromosome that will be 
%               included in the outut BED files. Must be 
%               maxChrNum >= minChrNum.
%               Default value: 22
%
% outdir
%               A string that specifies the path of the directory in which  
%               the output BED files are written. 
%               Default value "$INTERMEDIATE"
%
% ESIflag
%               Flag that determines whether this function performs 
%               computation of the entropic sensitivity index (ESI). 
%               0: no ESI computation. 
%               1: allow ESI computation.
%               Default value: 0
%
% MSIflag
%               Flag that determines whether this function performs 
%               computation of the methylation sensitivity index (MSI). 
%               0: no MSI computation. 
%               1: allow MSI computation.
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
%               The ratio regionSize/subregionSize must be an integer.
%               Default value: 150
%
% minNumCpG     The minimum number of CpG sites within an analysis 
%               subregion required for performing full methylation-based 
%               differential classification. 
%               Default value: 2
%
% thresh        
%
%               A scalar used as a threshold in methylation-based
%               differential classification.
%               Default value: 0.55
%
% threshDMU
%               A 1x6 vector containing threshold values used for
%               methylation-based differential classification.
%               Default value: [-1,-0.55,-0.1,0.1,0.55,1]
%
% threshDEU
%               A 1x8 vector containing threshold values used for
%               entropy-based differential classification.
%               Default value: [-1,-0.5,-0.3,-0.05,0.05,0.3,0.5,1]
%
% The default values of regionSize, subregionSize, minNumCpG, thresh, 
% threshDMU, and threshDEU should only be changed by an expert with 
% a detailed understanding of the code and the methods used. 
%

function makeBEDsForDiffMethAnalysis(prefix_1,prefix_2,analysis_path_1,analysis_path_2,...
	reference_path,varargin)

% Argument parser 
p = inputParser;

% Mandatory arguments
addRequired(p,'prefix_1')
addRequired(p,'prefix_2')
addRequired(p,'analysis_path_1')
addRequired(p,'analysis_path_2')
addRequired(p,'reference_path')

% Optional arguments
addParameter(p,'outdir',['..' filesep 'makeBedsForDiffMethAnalysis_out' filesep],...
               	@(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'minChrNum',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'maxChrNum',22,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'ESIflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MSIflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'MCflag',0,@(x)validateattributes(x,{'numeric'},{'nonempty','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))     
addParameter(p,'minNumCpG',2,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','nonnegative','scalar'})) 
addParameter(p,'thresh',0.55,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','positive','scalar'}))    
addParameter(p,'threshDMU',[-1,-0.55,-0.1,0.1,0.55,1],@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','nonnegative','scalar','size',[1,6]}))        
addParameter(p,'threshDEU',[-1,-0.5,-0.3,-0.05,0.05,0.3,0.5,1],@(x)validateattributes(x,...
		{'numeric'},{'nonempty','integer','nonnegative','scalar','size',[1,8]}))

% Parse optional arguments
parse(p,prefix_1,prefix_2,analysis_path_1,analysis_path_2,reference_path,varargin{:})
minChrNum       = p.Results.minChrNum;
maxChrNum       = p.Results.maxChrNum;
outdir          = p.Results.outdir;
ESIflag         = p.Results.ESIflag;
MSIflag         = p.Results.MSIflag;
MCflag          = p.Results.MCflag;
regionSize      = p.Results.regionSize;
subregionSize   = p.Results.subregionSize;
minNumCpG       = p.Results.minNumCpG;
thresh          = p.Results.thresh;
threshDMU       = p.Results.threshDMU;
threshDEU       = p.Results.threshDEU;

% Manual checks/corrections of inputs.
if analysis_path_1(end)~=filesep
    analysis_path_1=[analysis_path_1 filesep];
end
if analysis_path_2(end)~=filesep
    analysis_path_2=[analysis_path_2 filesep];
end
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if outdir(end)~=filesep
    outdir=[outdir filesep];
end

% Output BED format
formatSpec = '%s\t%u\t%u\t%f\n';

% Define hard coded options.

% dMML
hardCodedOptionsdMML = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dMML-';

% DMU
hardCodedOptionsDMU = 'track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-3.0:3.0 name=DMU-';

% dNME
hardCodedOptionsdNME = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dNME-';

% DEU
hardCodedOptionsDEU = 'track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-3.0:3.0 name=DEU-';

% JSD
hardCodedOptionsJSD = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=JSD-';

if ESIflag
	% dESI
	hardCodedOptionsdESI = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dESI-';
end

if MSIflag
	% dMSI
	hardCodedOptionsdMSI = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dMSI-';
end

if MCflag
	% dCAP
	hardCodedOptionsdCAP = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dCAP-';
	% dRDE
	hardCodedOptionsdRDE = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dRDE-';
end

% Check min and max chromosomes make sense
if minChrNum > maxChrNum
    fprintf(2,'error: minChrNum must be <= maxChrNum\n');
    return;
end

% Open BED files to write output.
% Proceed to open file handles.

% dMML
bedFileNamedMML = [outdir 'dMML-' prefix_1 '-VS-' prefix_2 '.bed'];
bedFileIDdMML   = fopen(bedFileNamedMML,'w'); 

% DMU
bedFileNameDMU  = [outdir 'DMU-' prefix_1 '-VS-' prefix_2 '.bed'];
bedFileIDDMU    = fopen(bedFileNameDMU,'w'); 

% dNME
bedFileNamedNME = [outdir 'dNME-' prefix_1 '-VS-' prefix_2 '.bed'];
bedFileIDdNME   = fopen(bedFileNamedNME,'w'); 

% DEU
bedFileNameDEU  = [outdir 'DEU-' prefix_1 '-VS-' prefix_2 '.bed'];
bedFileIDDEU    = fopen(bedFileNameDEU,'w'); 

% JSD
bedFileNameJSD  = [outdir 'JSD-' prefix_1 '-VS-' prefix_2 '.bed'];
bedFileIDJSD    = fopen(bedFileNameJSD,'w');

if ESIflag
    % dESI
	bedFileNamedESI = [outdir 'dESI-' prefix_1 '-VS-' prefix_2 '.bed'];
	bedFileIDdESI   = fopen(bedFileNamedESI,'w');
end

if MSIflag
    % dMSI
	bedFileNamedMSI = [outdir 'dMSI-' prefix_1 '-VS-' prefix_2 '.bed'];
	bedFileIDdMSI   = fopen(bedFileNamedMSI,'w');
end

if MCflag
	% dCAP
	bedFileNamedCAP = [outdir 'dCAP-' prefix_1 '-VS-' prefix_2 '.bed'];
	bedFileIDdCAP   = fopen(bedFileNamedCAP,'w'); 

	% dRDE
	bedFileNamedRDE = [outdir 'dRDE-' prefix_1 '-VS-' prefix_2 '.bed'];
	bedFileIDdRDE   = fopen(bedFileNamedRDE,'w'); 
end

% Write header information.

% dMML
fprintf(bedFileIDdMML,'%s%s-VS-%s\n',hardCodedOptionsdMML,prefix_1,prefix_2);

% DMU
fprintf(bedFileIDDMU,'%s%s-VS-%s\n',hardCodedOptionsDMU,prefix_1,prefix_2);

% dNME
fprintf(bedFileIDdNME,'%s%s-VS-%s\n',hardCodedOptionsdNME,prefix_1,prefix_2);

% DEU
fprintf(bedFileIDDEU,'%s%s-VS-%s\n',hardCodedOptionsDEU,prefix_1,prefix_2);

% JSD
fprintf(bedFileIDJSD,'%s%s-VS-%s\n',hardCodedOptionsJSD,prefix_1,prefix_2);

if ESIflag
	% dESI
	fprintf(bedFileIDdESI,'%s%s-VS-%s\n',hardCodedOptionsdESI,prefix_1,prefix_2);
end

if MSIflag
	% dMSI
	fprintf(bedFileIDdMSI,'%s%s-VS-%s\n',hardCodedOptionsdMSI,prefix_1,prefix_2);
end

if MCflag
	% dCAP
	fprintf(bedFileIDdCAP,'%s%s-VS-%s\n',hardCodedOptionsdCAP,prefix_1,prefix_2);
	% dRDE
	fprintf(bedFileIDdRDE,'%s%s-VS-%s\n',hardCodedOptionsdRDE,prefix_1,prefix_2); 
end

% Loop over chromosomes.

for chrNum = minChrNum:maxChrNum
    try
        chr_num_str = num2str(chrNum);
        chr_str     = ['chr' chr_num_str];

        % Find last CpG site on chromosome (to determine how many 
        % regions to loop through).       
        CpGdata = [reference_path 'CpGlocationChr' chr_num_str '.mat'];

        clear finalCpGloc CpGlocation          
        load(CpGdata,'finalCpGloc','CpGlocation');

        % Find all start base pairs for regions to be modeled 
        % on this chromosome.
        
        allStartBPs = int64(1):regionSize:int64(finalCpGloc); 
                        % Only need to go to the last CpG site, not to last
                        % base pair.

        % Try loading first file.
        clear mapObjData;
        analysis_file_1=[analysis_path_1 chr_str filesep prefix_1 '_analysis.mat'];
        try
            load(analysis_file_1,'mapObjData');
        	mapObjData1 = mapObjData;
        catch
            fprintf(2,['There was a problem loading file: ' analysis_file_1 ...
                '. Chromosome ' chr_num_str ' will not be included.\n']);
            continue
        end

        analysis_file_2=[analysis_path_2 chr_str filesep prefix_2 '_analysis.mat'];
        % Try loading second file.
        try
            load(analysis_file_2,'mapObjData');
        	mapObjData2 = mapObjData;
        catch
            fprintf(2,['There was a problem loading file: ' analysis_file_2 ...
                '. Chromosome ' chr_num_str ' will not be included.\n']);
            continue
        end
        clear mapObjData;
        
        % Loop through all regions on chromosome.
        for regionNum = 1:length(allStartBPs)
            try

                % Find key of this region in the hashtable.
                startBP = allStartBPs(regionNum);
                endBP   = startBP+regionSize-int64(1);
                locationPathName = [chr_str '/bp' num2str(startBP) ...
                                                     '-' num2str(endBP)];
                
                
                if isKey(mapObjData1,locationPathName) && ...
                                      isKey(mapObjData2,locationPathName)
                    
                    % Load structures for region (phenotype 1).
                    regionStruct1 = mapObjData1(locationPathName);
                    fullLProbs1   = regionStruct1.fullLProbs;
                    isModeled1    = regionStruct1.isModeled;
                    Ncg1          = regionStruct1.Ncg;
                    NME1          = regionStruct1.NME;
                   
                    if ESIflag
                        ESI1 = regionStruct1.ESI;
                        if isempty(ESI1)
                            fprintf(2,'Error: ESIflag = 1, but informME_run was run with ESIflag = 0.\n');
                            fprintf(2,'Error: In order to obtain ESI delete file:\n');
                            fprintf(2,analysis_file_1);
                            fprintf(2,'Error: and run informME_run again with --ESI flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdESI);
                            if MSIflag
                                fclose(bedFileIDdMSI);
                            end
                            if MCflag
                                fclose(bedFileIDdCAP);
                                fclose(bedFileIDdRDE);
                            end
                            return;
                        end
                    end

                    if MSIflag
                        MSI1 = regionStruct1.MSI;
                        if isempty(MSI1)
                            fprintf(2,'Error: MSIflag = 1, but informME_run was run with MSIflag = 0.\n');
                            fprintf(2,'Error: In order to obtain MSI delete file:\n');
                            fprintf(2,analysis_file_1);
                            fprintf(2,'Error: and run informME_run again with --MSI flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdMSI);
                            if ESIflag
                                fclose(bedFileIDdESI);
                            end
                            if MCflag
                                fclose(bedFileIDdCAP);
                                fclose(bedFileIDdRDE);
                            end
                            return;
                        end
                    end
                    
                    if MCflag
                        CAP1 = regionStruct1.CAP;
                        RDE1 = regionStruct1.RDE;
                        if isempty(CAP1) || isempty(RDE1)
                            fprintf(2,'Error: MCflag = 1, but informME_run was run with MCflag = 0.\n');
                            fprintf(2,'Error: In order to obtain MC delete file:\n');
                            fprintf(2,analysis_file_1);
                            fprintf(2,'Error: and run informME_run again with --MC flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdCAP);
                            fclose(bedFileIDdRDE);
                            if ESIflag
                                fclose(bedFileIDdESI);
                            end
                            if MSIflag
                                fclose(bedFileIDdMSI);
                            end
                            return;
                        end
                    end
                    
                    % Load structures for region (penotype 2).
                    regionStruct2 = mapObjData2(locationPathName);
                    fullLProbs2   = regionStruct2.fullLProbs;
                    isModeled2    = regionStruct2.isModeled;
                    Ncg2          = regionStruct2.Ncg;
                    NME2          = regionStruct2.NME;

                    if ESIflag
                        ESI2 = regionStruct2.ESI;
                        if isempty(ESI2)
                            fprintf(2,'Error: ESIflag = 1, but informME_run was run with ESIflag = 0.\n');
                            fprintf(2,'Error: In order to obtain ESI delete file:\n');
                            fprintf(2,analysis_file_2);
                            fprintf(2,'Error: and run informME_run again with --ESI flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdESI);
                            if MSIflag
                                fclose(bedFileIDdMSI);
                            end
                            if MCflag
                                fclose(bedFileIDdCAP);
                                fclose(bedFileIDdRDE);
                            end
                            return;
                        end
                    end

                    if MSIflag
                        MSI2 = regionStruct2.MSI;
                        if isempty(MSI2)
                            fprintf(2,'Error: MSIflag = 1, but informME_run was run with MSIflag = 0.\n');
                            fprintf(2,'Error: In order to obtain MSI delete file:\n');
                            fprintf(2,analysis_file_2);
                            fprintf(2,'Error: and run informME_run again with --MSI flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdMSI);
                            if ESIflag
                                fclose(bedFileIDdESI);
                            end
                            if MCflag
                                fclose(bedFileIDdCAP);
                                fclose(bedFileIDdRDE);
                            end
                            return;
                        end
                    end
                    
                    if MCflag
                        CAP2 = regionStruct2.CAP;
                        RDE2 = regionStruct2.RDE;
                        if isempty(CAP2) || isempty(RDE2)
                            fprintf(2,'Error: MCflag = 1, but informME_run was run with MCflag = 0.\n');
                            fprintf(2,'Error: In order to obtain MC delete file:\n');
                            fprintf(2,analysis_file_2);
                            fprintf(2,'Error: and run informME_run again with --MC flag.\n');
                            fclose(bedFileIDdMML);
                            fclose(bedFileIDDMU);
                            fclose(bedFileIDdNME);
                            fclose(bedFileIDDEU);
                            fclose(bedFileIDJSD);
                            fclose(bedFileIDdCAP);
                            fclose(bedFileIDdRDE);
                            if ESIflag
                                fclose(bedFileIDdESI);
                            end
                            if MSIflag
                                fclose(bedFileIDdMSI);
                            end
                            return;
                        end
                    end

                    if sum(Ncg1~=Ncg2)
                        fprintf(2,'ERROR: number of CpG sites are different in the two phenotypes\n');
                        ME = MException('VerifyOutput:OutOfBounds','Results are outside allowable limits');
			fclose(bedFileIDdMML);
                        fclose(bedFileIDDMU);
			fclose(bedFileIDdNME);
			fclose(bedFileIDDEU);
			fclose(bedFileIDJSD);
                        if ESIflag
                            fclose(bedFileIDdESI);
                        end
                        if MSIflag
                            fclose(bedFileIDdMSI);
                        end
                        if MCflag
                            fclose(bedFileIDdCAP);
                            fclose(bedFileIDdRDE);
                        end
                        throw(ME);
                    end
                    
                    % Compute regional quantities.
                    [lower_index,upper_index] =...
                             findSortedIndices(CpGlocation,startBP,endBP);
                    CpGlocInReg = CpGlocation(lower_index:upper_index);
                    numCpGinReg = length(CpGlocInReg);
                   
                    subRegCount = 0;
                    for offset = 0:subregionSize:(regionSize-1)
                        subRegCount = subRegCount + 1;
                        subStartBP  = startBP + offset;
                        subEndBP    = min(subStartBP + subregionSize...
                                                - int64(1),finalCpGloc);
                        
                        numCpGsInSubReg = Ncg1(subRegCount);
                        
                        if numCpGsInSubReg>0 && (isModeled1(subRegCount)) ...
                                && (isModeled2(subRegCount)) 
                                 % Subregion is modeled in both phenotypes. 
                            
                            [lower_index,upper_index] = ...
                                findSortedIndices(CpGlocInReg,subStartBP,subEndBP);
                            if numCpGsInSubReg==1 && ...
                                    (lower_index==1||upper_index==numCpGinReg) 
                                % Only 1 CpG site and that is at boundary 
                                % of a genomic region used for parameter 
                                % estimation.
                                continue; % Move onto next subregion.
                            end

                            % Print dNME, dESI, dMSI, dCAP, dRDE
                            
                            dNMERegion = NME1(subRegCount)-NME2(subRegCount);
                            
                            if dNMERegion < inf
                                fprintf(bedFileIDdNME,formatSpec,...
                                              chr_str,int64(subStartBP-1),...
                                              int64(subEndBP-1),dNMERegion);
                            end
                            
                            if ESIflag
                                dESI = ESI1(subRegCount)-ESI2(subRegCount);
                                if dESI < inf
                                    fprintf(bedFileIDdESI,formatSpec,...
                                        chr_str,int64(subStartBP-1),...
                                        int64(subEndBP-1),dESI);
                                end
                            end

                            if MSIflag
                                dMSI = MSI1(subRegCount)-MSI2(subRegCount);
                                if dMSI < inf
                                    fprintf(bedFileIDdMSI,formatSpec,...
                                        chr_str,int64(subStartBP-1),...
                                        int64(subEndBP-1),dMSI);
                                end
                            end
                                   
                            if MCflag
                                dCAP = CAP1(subRegCount)-CAP2(subRegCount);
                                if dCAP < inf
                                    fprintf(bedFileIDdCAP,formatSpec,...
                                        chr_str,int64(subStartBP-1),...
                                        int64(subEndBP-1),dCAP);
                                end
                                
                                dRDE = RDE1(subRegCount)-RDE2(subRegCount);
                                if dRDE < inf
                                    fprintf(bedFileIDdRDE,formatSpec,...
                                        chr_str,int64(subStartBP-1),...
                                        int64(subEndBP-1),dRDE);
                                end
                            end
	            
                            % Perform entropy-based differential classification. 
                            if dNMERegion >= threshDEU(7) 
                                                % Strongly hyperentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),3);
                            elseif dNMERegion >= threshDEU(6) 
                                                % Moderately hyperentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),2);
                            elseif dNMERegion >= threshDEU(5) 
                                                % Weakly hyperentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),1);
                            elseif dNMERegion > threshDEU(4) 
                                                % Isoentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),0);
                            elseif dNMERegion > threshDEU(3) 
                                                 % Weakly hypoentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),-1);
                            elseif dNMERegion > threshDEU(2) 
                                                 % Moderately hypoentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),-2);
                            elseif dNMERegion >=threshDEU(1) 
                                                 % Strongly hypoentropic.
                                fprintf(bedFileIDDEU,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),-3);
                            end % End entropy-based differential classification. 
                            
                            % Compute probability distribution of the 
                            % difference D = L1-L2 in methylation levels 
                            % between the two phenotypes as well as the 
                            % mean of this distribution (which is the dMML).
                            pL1 = full(fullLProbs1(subRegCount,1:(numCpGsInSubReg+1)));
                            pL2 = full(fullLProbs2(subRegCount,1:(numCpGsInSubReg+1)));
                            
                            pD    = conv(pL1,pL2(end:-1:1));
                            Dvals = -1:(1/double(numCpGsInSubReg)):1;
                            
                            % Compute and print dMML.
                            dMML = dot(pD,Dvals);
                            if dMML < inf
                                fprintf(bedFileIDdMML,formatSpec,...
                                              chr_str,int64(subStartBP-1),...
                                              int64(subEndBP-1),dMML);
                            end
                                                 
                            % Compute and print JSD.
                            
                            pL1 = pL1(:); % Ensures column vector.
                            pL2 = pL2(:); % Ensures column vector.
                           
                            JSD = sqrt(sum(pL1(pL1>0).*log2(2.*pL1(pL1>0)...
                                           ./(pL1(pL1>0)+pL2(pL1>0))))/2 ...
                                     + sum(pL2(pL2>0).*log2(2.*pL2(pL2>0)...
                                           ./(pL1(pL2>0)+pL2(pL2>0))))/2);
                            
                            if JSD < inf      
                                fprintf(bedFileIDJSD,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),JSD);
                            end
                                                                           
                            % Compute coarse probabilities q of D.
                            
                            qProbsRegion    = zeros(5,1);
                            
                            qProbsRegion(1) = sum(pD((threshDMU(1)<=Dvals)...
                                                     &(Dvals<=threshDMU(2)))); 
                            qProbsRegion(2) = sum(pD((threshDMU(2)<Dvals)...
                                                     &(Dvals<=threshDMU(3))));  
                            qProbsRegion(3) = sum(pD((threshDMU(3)<Dvals)...
                                                      &(Dvals<threshDMU(4))));   
                            qProbsRegion(4) = sum(pD((threshDMU(4)<=Dvals)...
                                                      &(Dvals<threshDMU(5))));  
                            qProbsRegion(5) = sum(pD((threshDMU(5)<=Dvals)...
                                                     &(Dvals<=threshDMU(6)))); 
                                                 
                            % Perform methylation-based differential classification.
                            
                            if numCpGsInSubReg >= minNumCpG % For more than 
                                                            % minNumCpG sites. 
                                                 
                                if qProbsRegion(3) > thresh
                                % Phenotypes 1 & 2 are isomethylated.
                                    fprintf(bedFileIDDMU,formatSpec,...
                                                 chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),0);
                                                    
                                elseif ((qProbsRegion(1)+qProbsRegion(2))...
                                               /(1-qProbsRegion(3))) > thresh
                                           
                                    if qProbsRegion(1) > thresh
                                    % Phenotype 1 is strongly hypomethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-3);
                                                   
                                    elseif (qProbsRegion(1)+qProbsRegion(2)) > thresh
                                    % Phenotype 1 is moderately hypomethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                     chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-2);
                                                   
                                    else 
                                    % Phenotype 1 is weakly hypomethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                     chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-1);
				
                                    end

                                elseif ((qProbsRegion(4)+qProbsRegion(5)) ...
                                               /(1-qProbsRegion(3))) > thresh
                                           
                                    if qProbsRegion(5) > thresh
                                    % Phenotype 1 is strongly hypermethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                     chr_str,int64(subStartBP-1),...
                                                           int64(subEndBP-1),3);
                                                       
                                    elseif (qProbsRegion(4)+qProbsRegion(5)) > thresh
                                    % Phenotype 1 is moderately hypermethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                     chr_str,int64(subStartBP-1),...
                                                           int64(subEndBP-1),2);
                                                       
                                    else 
                                    % Phenotype 1 is weakly hypermethylated.
                                        fprintf(bedFileIDDMU,formatSpec,...
                                                     chr_str,int64(subStartBP-1),...
                                                           int64(subEndBP-1),1);
                                    end
                                    
                                end

                            else % For less than minNumCpG sites.
                                if qProbsRegion(3) > thresh
                                    
                                    % Phenotypes 1 & 2 are isomethylated.
                                    fprintf(bedFileIDDMU,formatSpec,...
                                                 chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),0);
                                                    
                                elseif (qProbsRegion(4)+qProbsRegion(5)) > thresh
                                    % Phenotype 1 is moderately hypermethylated.
                                    fprintf(bedFileIDDMU,formatSpec,...
                                                 chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),2);
                                                    
                                elseif (qProbsRegion(1)+qProbsRegion(2)) > thresh
                                     % Phenotype 1 is moderately hypomethylated.
                                     fprintf(bedFileIDDMU,formatSpec,...
                                                  chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),-2);
                                                    
                                elseif (qProbsRegion(4)+qProbsRegion(5))...
                                               /(1-qProbsRegion(3)) > thresh
                                    % Phenotype 1 is weakly hypermethylated.
                                    fprintf(bedFileIDDMU,formatSpec,...
                                                 chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),1);
                                                    
                                 elseif (qProbsRegion(1)+qProbsRegion(2))...
                                               /(1-qProbsRegion(3)) > thresh 
                                    % Phenotype 1 is weakly hypomethylated.
                                    fprintf(bedFileIDDMU,formatSpec,...
                                                 chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-1);
                                                   
                                end
                                    
                            end % End methylation-based differential classification.
                          
                        end % End subregion modeled. 
                         
                    end % End loop over subregions.
                    
                end % Data in both regions.
                
            catch ME1
                fprintf(2,['Error at Chr' chr_num_str ':bp' num2str(startBP) ...
                                        '-' num2str(endBP) ' with message:\n%s'],ME1.message)
            end
        end % End loop over regions.
    catch ME2
        fprintf(2,['Error at Chr ' num2str(chrNum) ' with message:\n%s'],ME2.message)
    end
    
end % End loop over chromosomes.

% Close all files currently open.
fclose(bedFileIDdMML);
fclose(bedFileIDDMU);
fclose(bedFileIDdNME);
fclose(bedFileIDDEU);
fclose(bedFileIDJSD);
if ESIflag
	fclose(bedFileIDdESI);
end
if MSIflag
	fclose(bedFileIDdMSI);
end
if MCflag
	fclose(bedFileIDdCAP);
	fclose(bedFileIDdRDE);
end

