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
%%%%%%%%                makeBEDsForMethAnalysis.m                  %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%                Last Modified: 12/08/2016                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function makes BED files for the methylation analysis results 
% obtained by means of MethAnalysisForChr.m for a single phenotype. 
%
% USAGE (default):
%
% makeBedsForMethAnalysis(prefix,analysis_path,reference_path)
%
% USAGE (optional):
%
% Example of optional usage with additional input parameters.
% makeBedsForMethAnalysis(prefix,analysis_path,reference_path,'outdir',
%	'/path/to/output')
% 
% MANDATORY INPUTS:
%
% prefix
%               A string that specifies the name of the phenotype. 
%
% analysis_path
%               A string that specifies the path of the directory in which  
%               the model was constructed. 
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
%               included in the output BED files. Must be 
%               maxChrNum >= minChrNum.
%               Default value: 22
%
% outdir
%               A string that specifies the path of the directory in which  
%               the output BED files are written. 
%               Default value './'
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
% thresh
%               A scalar used as a threshold in methylation-based
%               classification.
%               Default value: 0.4
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
% The default values of thresh, regionSize, and subregionSize should only  
% be changed by an expert with a detailed understanding of the code and 
% the methods used. 
%

function makeBedsForMethAnalysis(prefix,analysis_path,reference_path,varargin)

% Argument parser
p = inputParser;

% Mandatory arguments
addRequired(p,'prefix')
addRequired(p,'analysis_path')
addRequired(p,'reference_path')

% OPtional arguments
addParameter(p,'outdir',['.' filesep 'singleMethAnalysisToBed_out' filesep],...
               	@(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'minChrNum',1,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'maxChrNum',22,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'ESIflag',0,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','scalar'}))
addParameter(p,'MCflag',0,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','scalar'}))
addParameter(p,'thresh',0.4,@(x)validateattributes(x,{'numeric'},...
		{'nonempty','positive','scalar'}))
addParameter(p,'regionSize',int64(3000),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),@(x)validateattributes(x,{'numeric'},...
		{'nonempty','integer','positive','scalar'}))

% Parse optional arguments
parse(p,prefix,analysis_path,reference_path,varargin{:})
minChrNum       = p.Results.minChrNum;
maxChrNum       = p.Results.maxChrNum;
outdir          = p.Results.outdir;
ESIflag         = p.Results.ESIflag;
MCflag          = p.Results.MCflag;
thresh          = p.Results.thresh;
regionSize      = p.Results.regionSize;
subregionSize   = p.Results.subregionSize;

% Manual checks/corrections of inputs.
if analysis_path(end)~=filesep
    analysis_path=[analysis_path filesep];
end
if reference_path(end)~=filesep
    reference_path=[reference_path filesep];
end
if outdir(end)~=filesep
    outdir=[outdir filesep];
end

% Define hard coded options.
% MML
hardCodedOptionsMML = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=MML-';
% NME
hardCodedOptionsNME = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=NME-';
% METH
hardCodedOptionsMETH = 'track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-2.0:2.0 name=METH-';
% VAR
hardCodedOptionsVAR = 'track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=0.0:3.0 name=VAR-';
% ENTR
hardCodedOptionsENTR = 'track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-2.0:2.0 name=ENTR-';

if ESIflag
	% ESI
	hardCodedOptionsESI = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=ESI-';
end

if MCflag
    	% TURN
	hardCodedOptionsTURN = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=TURN-';
    	% CAP
	hardCodedOptionsCAP = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=CAP-';
	% RDE
	hardCodedOptionsRDE = 'track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=RDE-';
end

% Check min and max chromosomes make sense
if minChrNum > maxChrNum
    fprintf(2,'error: minChrNum must be <= maxChrNum\n');
    return;
end

% Open BED files to write output.
% Proceed to open file handles
% MML
bedFileNameMML = [outdir 'MML-' prefix '.bed'];
bedFileIDMML    = fopen(bedFileNameMML,'w'); 
% NME
bedFileNameNME = [outdir 'NME-' prefix '.bed'];
bedFileIDNME   = fopen(bedFileNameNME,'w'); 
% METH
bedFileNameMETH = [outdir 'METH-' prefix '.bed'];
bedFileIDMETH   = fopen(bedFileNameMETH,'w'); 
% VAR
bedFileNameVAR = [outdir 'VAR-' prefix '.bed'];
bedFileIDVAR   = fopen(bedFileNameVAR,'w'); 
% ENTR
bedFileNameENTR = [outdir 'ENTR-' prefix '.bed'];
bedFileIDENTR   = fopen(bedFileNameENTR,'w'); 
% ESI
if ESIflag
	bedFileNameESI = [outdir 'ESI-' prefix '.bed'];
	bedFileIDESI   = fopen(bedFileNameESI,'w'); 
end

% MC
if MCflag   
    % TURN
	bedFileNameTURN = [outdir 'TURN-' prefix '.bed'];
	bedFileIDTURN   = fopen(bedFileNameTURN,'w'); 
	% CAP
	bedFileNameCAP = [outdir 'CAP-' prefix '.bed'];
	bedFileIDCAP   = fopen(bedFileNameCAP,'w'); 
	% RDE
	bedFileNameRDE = [outdir 'RDE-' prefix '.bed'];
	bedFileIDRDE   = fopen(bedFileNameRDE,'w');
end

% Write header information.
% MML
fprintf(bedFileIDMML,'%s%s\n',hardCodedOptionsMML,prefix);
% NME
fprintf(bedFileIDNME,'%s%s\n',hardCodedOptionsNME,prefix);
% METH
fprintf(bedFileIDMETH,'%s%s\n',hardCodedOptionsMETH,prefix);
% VAR
fprintf(bedFileIDVAR,'%s%s\n',hardCodedOptionsVAR,prefix);
% ENTR
fprintf(bedFileIDENTR,'%s%s\n',hardCodedOptionsENTR,prefix);

% ESI
if ESIflag
	fprintf(bedFileIDESI,'%s%s\n',hardCodedOptionsESI,prefix);    
end

% MC
if MCflag
  	% TURN
	fprintf(bedFileIDTURN,'%s%s\n',hardCodedOptionsTURN,prefix);
	% CAP
	fprintf(bedFileIDCAP,'%s%s\n',hardCodedOptionsCAP,prefix);
	% RDE
	fprintf(bedFileIDRDE,'%s%s\n',hardCodedOptionsRDE,prefix);
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

        % Find all start base pairs for regions to be modeled on this chromosome.
        % Only need to go to the last CpG site, not to last base pair.
        allStartBPs = int64(1):regionSize:int64(finalCpGloc); 

        % Load relevant result files
        clear mapObjData;
		analysis_file=[analysis_path chr_str filesep prefix '_analysis.mat'];
		try
        	load(analysis_file,'mapObjData');
		catch
			fprintf(2,['There was a problem loading file: ' analysis_file ...
				'. Chromosome ' chr_num_str ' will not be included.\n']); 
			continue
		end
        
        % Loop through all regions on chromosome. 
        for regionNum = 1:length(allStartBPs)
            try
                
                % Find key of this region in the hashtable.
                startBP = allStartBPs(regionNum);
                endBP   = startBP+regionSize-int64(1);
                locationPathName = [chr_str '/bp' num2str(startBP) '-' num2str(endBP)];
                
                if isKey(mapObjData,locationPathName)
                    % Load structures for region.
                    regionStruct = mapObjData(locationPathName);
                    cLProbs      = regionStruct.cLProbs;  
                    isModeled    = regionStruct.isModeled;
                    Ncg          = regionStruct.Ncg;
                    NME          = regionStruct.NME;
                    MML          = regionStruct.MML;
                    
                    if ESIflag
                        ESI = regionStruct.ESI;
                        if isempty(ESI)
                            fprintf(2,'Error: ESIflag = 1, but the previous methylation analysis step was run with ESIflag = 0.\n');
                            fprintf(2,'Rerun this MATLAB function with ESIflag = 0 or rerun the previous methylation analysis step with ESIflag = 1.\n');
                            fclose(bedFileIDMML);
                            fclose(bedFileIDNME);
                            fclose(bedFileIDMETH);
                            fclose(bedFileIDVAR);
                            fclose(bedFileIDENTR);
                            if ESIflag
                                fclose(bedFileIDESI);
                            end
                            if MCflag
                                fclose(bedFileIDTURN);
                                fclose(bedFileIDCAP);
                                fclose(bedFileIDRDE);
                            end
                            return;
                        end
                    end
                    
                    if MCflag
                        TURN = regionStruct.TURN;
                        CAP  = regionStruct.CAP;
                        RDE  = regionStruct.RDE;
                        if isempty(CAP) || isempty(RDE) || isempty(TURN)
                            fprintf(2,'Error: MCflag = 1, but the previous methylation analysis step was run with MCflag = 0.\n');
                            fprintf(2,'Rerun this MATLAB function with MCflag = 0 or rerun the previous methylationm analysis step with MCflag = 1.\n');
                            fclose(bedFileIDMML);
                            fclose(bedFileIDNME);
                            fclose(bedFileIDMETH);
                            fclose(bedFileIDVAR);
                            fclose(bedFileIDENTR);
                            if ESIflag
                                fclose(bedFileIDESI);
                            end
                            if MCflag
                                fclose(bedFileIDTURN);
                                fclose(bedFileIDCAP);
                                fclose(bedFileIDRDE);
                            end
                            return;
                        end
                    end
                   
                    % Compute regional quantities.                
                    [lower_index,upper_index] = findSortedIndices(CpGlocation,startBP,endBP);
                    CpGlocInReg = CpGlocation(lower_index:upper_index);
                    numCpGinReg = length(CpGlocInReg);
                    
                    subRegCount = 0;
                    for offset = 0:subregionSize:(regionSize-1)
                        subRegCount = subRegCount + 1;
                        subStartBP  = startBP + offset;
                        subEndBP    = min(subStartBP + subregionSize - int64(1),finalCpGloc);
                     
                        numCpGsInSubReg = Ncg(subRegCount);
                        
                        if numCpGsInSubReg>0 && isModeled(subRegCount) 
                                                       % Subregion is modeled. 
                            
                            [lower_index,upper_index] = findSortedIndices(CpGlocInReg,subStartBP,subEndBP);
                            if numCpGsInSubReg==1 && (lower_index==1||upper_index==numCpGinReg) 
                                % Only 1 CpG site and that is at boundary
                                % of a genomic region used for parameter 
                                % estimation. 
                                continue; % Move onto next subregion.
                            end
                            
                            cLProbsRegion = cLProbs(subRegCount,:);
                            MMLregion     = MML(subRegCount);
                            NMEregion     = NME(subRegCount);
                            
                            if ESIflag
                                ESIregion = ESI(subRegCount);
                            end
                            
                            if MCflag
                                TURNregion = TURN(subRegCount);
                                CAPregion  = CAP(subRegCount);
                                RDEregion  = RDE(subRegCount); 
                            end
                            
                            % Print MML, ESI, TURN, CAP, and RDE.
			    formatSpec = '%s\t%u\t%u\t%f\n';
                            if MMLregion<inf
                                fprintf(bedFileIDMML,formatSpec,...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),MMLregion);
                            end
                            
                            if  ESIflag && ESIregion<inf
                                fprintf(bedFileIDESI,'%s %12u %12u %f\n',...
                                              chr_str,int64(subStartBP-1),...
                                              int64(subEndBP-1),ESIregion);
                            end
                            
                            if MCflag && TURNregion<inf
                                fprintf(bedFileIDTURN,'%s %12u %12u %f\n',...
                                              chr_str,int64(subStartBP-1),...
                                              int64(subEndBP-1),TURNregion);
                            end
                            
                            if MCflag && CAPregion<inf
                                fprintf(bedFileIDCAP,'%s %12u %12u %f\n',...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),CAPregion);
                            end
                            
                            if MCflag && RDEregion<inf
                                fprintf(bedFileIDRDE,'%s %12u %12u %f\n',...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),RDEregion);
                            end
                            
                            % Print NME and perform NME classification.
                            if NMEregion<inf
                                fprintf(bedFileIDNME,'%s %12u %12u %f\n',...
                                             chr_str,int64(subStartBP-1),...
                                             int64(subEndBP-1),NMEregion);
                                if NMEregion >= 0.99 
                                    fprintf(bedFileIDENTR,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),2);
                                elseif NMEregion >= 0.92 
                                    fprintf(bedFileIDENTR,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),1);
                                elseif NMEregion > 0.44 
                                    fprintf(bedFileIDENTR,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),0);
                                elseif NMEregion > 0.28
                                    fprintf(bedFileIDENTR,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),-1);
                                elseif NMEregion >= 0    
                                    fprintf(bedFileIDENTR,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),-2);
                                end
                            end 
                            
                            % Perform MML classification.
                            if (numCpGsInSubReg > 1) % More than one CpG site
                                
                                unmethProb = cLProbsRegion(1)+cLProbsRegion(2);

                                if unmethProb > 1-thresh % Unmethylated.
                                    if cLProbsRegion(1) > 1-thresh
                                                  % Highly unmethylated.
                                        fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),-2);
                                    else          % Partially Unmethylated.
                                        fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),-1);
                                    end
                                elseif unmethProb < thresh % Methylated. 
                                    if cLProbsRegion(4) > 1-thresh
                                                    % Highly methylated.
                                        fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),2);
                                    else            % Partially Methylated.
                                        fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),1);
                                    end
                                elseif ~(sum(isnan(cLProbsRegion))>0) % Variable.
                                    fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),0); 
                                                       % Print 0 in METH track.
                                                       
                                    pratio1 = cLProbsRegion(1)/unmethProb;
                                    pratio2 = cLProbsRegion(4)/(1-unmethProb);
                                    if (pratio1 <= thresh) && (pratio2 <= thresh)
                                                             % Mixed.
                                        fprintf(bedFileIDVAR,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),1);
                                    elseif (pratio1 < 1-thresh) && (pratio2 < 1-thresh)
                                                            % Highly mixed.
                                        fprintf(bedFileIDVAR,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),2);
                                    elseif (pratio1 >= 1-thresh) && (pratio2 >= 1-thresh)
                                                             % Bistable.
                                        fprintf(bedFileIDVAR,'%s %12u %12u %f\n',...
                                            chr_str,int64(subStartBP-1),...
                                            int64(subEndBP-1),3);
                                    end % End variable classification.
                                end % End unmethylated/methylated/variable classification
                            else % Only one CpG site.
                                unmethProb = cLProbsRegion(1)+cLProbsRegion(2);
                                if unmethProb > 1-thresh    % Unmethylated.
                                              
                                    fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),-1);
                                elseif unmethProb < thresh   % Methylated.
                                    fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),1); 
                                elseif ~(sum(isnan(cLProbsRegion))>0) % Mixed. 
                                    fprintf(bedFileIDMETH,'%s %12u %12u %f\n',...
                                                  chr_str,int64(subStartBP-1),...
                                                  int64(subEndBP-1),0); 
                                    fprintf(bedFileIDVAR,'%s %12u %12u %f\n',...
                                                 chr_str,int64(subStartBP-1),...
                                                 int64(subEndBP-1),1);
                                end
                            end % End methylation-based classificationn.                           
                        end % End subregion modeled.
                    end % End loop over subregions.
                end % End data present in region.
            catch ME1
                fprintf(2,['Error at Chr' chr_num_str ':bp' num2str(startBP) ...
                                        '-' num2str(endBP) ' with message:\n%s\n'],ME1.message)
            end
        end % End loop over regions.
    catch ME2
        fprintf(2,['Error at Chr ' num2str(chrNum) ' with message: \n%s\n'],ME2.message)
    end
end % End loop over chromosomes.

% Close all files currently open.
fclose(bedFileIDMML);
fclose(bedFileIDNME);
fclose(bedFileIDMETH);
fclose(bedFileIDVAR);
fclose(bedFileIDENTR);

if ESIflag
	fclose(bedFileIDESI);
end

if MCflag
    fclose(bedFileIDTURN);
	fclose(bedFileIDCAP);
	fclose(bedFileIDRDE);
end

end
