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
%%%%%%%%%%%             Last Modified: 05/28/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function makes the relevant bed files for a given phenotype. It
% assumes that statistical estimation has been completed using
% EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
% with AnalysisForChr.m/MergeAnalysis.m for the phenotype.
%
% Example Default Usage:
% MakeBEDfilesAnalysis(phenoName)
%
% Example Usage Modifying optional parameter "maxChrNum" as name-value pair:
% MakeBEDfilesAnalysis(phenoName,'maxChrNum',5)
% 
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% phenoName
%               A string with the name of the 
%               phenotype that is has already been modeled by
%               EstParamsForChr.m
%
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MakeBEDfilesAnalysis(...,'inputName',inputValue):
%
% species
%               A string detailing the species from which the data is
%               obtained. i.e., 'Human' or 'Mouse'. Default: 'Human'
%
% minChrNum
%               A number specifying the starting chromosome to output to bed
%               file. So total chr to go out to file will be
%               minChrNum:maxChrNum. Default value: 1.
%
% maxChrNum
%               A number specifying how many chromosomes to output to bed
%               file. So total chr to go out to file will be
%               minChrNum:maxChrNum. Default value: 22.
%
% resultsPathRoot
%               A string with the path to the results directory where the
%               estimation and analysis results are written in its
%               subdirectories. Default value './results/'
%
% genomePathRoot
%               A string with the path to the results directory where the
%               genome analysis results are written in its subdirectories.
%               Default value: '../ParseBAMfile/'
%
% Any other optional inputs should only be modified by professionals with a
% detailed understanding of the code.

function MakeBEDfilesAnalysis(phenoName,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'phenoName')
addParameter(p,'species','Human',...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'minChrNum',1,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'maxChrNum',22,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'regionSize',int64(3000),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'subregionSize',int64(150),...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))
addParameter(p,'kappa',1.5,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','positive','scalar'}))

parse(p,phenoName,varargin{:})

species         = p.Results.species;
minChrNum       = p.Results.minChrNum;
maxChrNum       = p.Results.maxChrNum;
resultsPathRoot = p.Results.resultsPathRoot;
genomePathRoot  = p.Results.genomePathRoot;
regionSize      = p.Results.regionSize;
subregionSize   = p.Results.subregionSize;
kappa           = p.Results.kappa;

%
% Manual checks/corrections of inputs
%

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end
if resultsPathRoot(end)~=filesep
    resultsPathRoot=[resultsPathRoot filesep];
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hard coded options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%An
hardCodedOptionsAn = 'track type=bedGraph visibility=full autoScale=on name=An-';

%Cn
hardCodedOptionsCn = 'track type=bedGraph visibility=full autoScale=on name=Cn-';

%MNML
hardCodedOptionsMNML = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=MNL-';

%Monostable
hardCodedOptionsMonostable = 'track type=bedGraph visibility=full autoScale=off viewLimits=-2.0:2.0 name=METH-';

%Bistable
hardCodedOptionsBistable = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:3.0 name=MIXT-';

%NMEyClass
hardCodedOptionsNMEyClass = 'track type=bedGraph visibility=full autoScale=off viewLimits=-2.0:2.0 name=ENTR-';

%NMEy
hardCodedOptionsNMEy = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=NME-';

%MST
hardCodedOptionsMST = 'track type=bedGraph visibility=full autoScale=on alwaysZero=on name=MDS-';

%MML
hardCodedOptionsMML = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=MDL-';

%ESI
hardCodedOptionsESI = 'track type=bedGraph visibility=full autoScale=on alwaysZero=on name=ESI-';

%NCU
hardCodedOptionsNCU = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=CAP-';

%NRDE
hardCodedOptionsNRDE = 'track type=bedGraph visibility=full autoScale=on alwaysZero=on name=RDE-';

%PCC
%hardCodedOptionsPCC = 'track type=bedGraph visibility=full autoScale=off viewLimits=-1.0:1.0 name=PCC-';

%TURN
hardCodedOptionsTURN = 'track type=bedGraph visibility=full autoScale=on alwaysZero=on name=TURN-';

%MENT
hardCodedOptionsMENT = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=MENT-';

%XENT
hardCodedOptionsXENT = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=XENT-';

% 
% %ALPHA
% hardCodedOptionsALPHA = 'track type=bedGraph visibility=full autoScale=on name=ALPHA-';
% 
% %BETA
% hardCodedOptionsBETA = 'track type=bedGraph visibility=full autoScale=on name=BETA-';
% 
% %GAMMA
% hardCodedOptionsGAMMA = 'track type=bedGraph visibility=full autoScale=on name=GAMMA-';

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open .bed files to write output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%check that bedfile folder exists, create it if not
if ~exist([resultsPathRoot species filesep 'BEDFILES'],'dir')
    mkdir([resultsPathRoot species filesep 'BEDFILES']);
end

if minChrNum > maxChrNum
    disp('error: minChrNum must be <= maxChrNum');
    return;
end


%An
bedFileNameAn = [resultsPathRoot species filesep 'BEDFILES' filesep 'An-' phenoName '.bed'];
bedFileIDAn   = fopen(bedFileNameAn,'w'); 

%Monostable
bedFileNameCn = [resultsPathRoot species filesep 'BEDFILES' filesep 'Cn-' phenoName '.bed'];
bedFileIDCn   = fopen(bedFileNameCn,'w'); 

%Monostable
bedFileNameMNML = [resultsPathRoot species filesep 'BEDFILES' filesep 'MNL-' phenoName '.bed'];
bedFileIDMNML   = fopen(bedFileNameMNML,'w'); 

%Monostable
bedFileNameMonostable = [resultsPathRoot species filesep 'BEDFILES' filesep 'METH-' phenoName '.bed'];
bedFileIDMonostable   = fopen(bedFileNameMonostable,'w'); 

%Bistable
bedFileNameBistable = [resultsPathRoot species filesep 'BEDFILES' filesep 'MIXT-' phenoName '.bed'];
bedFileIDBistable   = fopen(bedFileNameBistable,'w'); 

%NMEyClass
bedFileNameNMEyClass = [resultsPathRoot species filesep 'BEDFILES' filesep 'ENTR-' phenoName '.bed'];
bedFileIDNMEyClass   = fopen(bedFileNameNMEyClass,'w'); 

%NMEy
bedFileNameNMEy = [resultsPathRoot species filesep 'BEDFILES' filesep 'NME-' phenoName '.bed'];
bedFileIDNMEy   = fopen(bedFileNameNMEy,'w'); 

%MST
bedFileNameMST = [resultsPathRoot species filesep 'BEDFILES' filesep 'MDS-' phenoName '.bed'];
bedFileIDMST   = fopen(bedFileNameMST,'w'); 

%MML
bedFileNameMML = [resultsPathRoot species filesep 'BEDFILES' filesep 'MDL-' phenoName '.bed'];
bedFileIDMML   = fopen(bedFileNameMML,'w'); 

%ESI
bedFileNameESI = [resultsPathRoot species filesep 'BEDFILES' filesep 'ESI-' phenoName '.bed'];
bedFileIDESI   = fopen(bedFileNameESI,'w'); 

%NCU
bedFileNameNCU = [resultsPathRoot species filesep 'BEDFILES' filesep 'CAP-' phenoName '.bed'];
bedFileIDNCU   = fopen(bedFileNameNCU,'w'); 

%NRDE
bedFileNameNRDE = [resultsPathRoot species filesep 'BEDFILES' filesep 'RDE-' phenoName '.bed'];
bedFileIDNRDE   = fopen(bedFileNameNRDE,'w'); 

%TURN
bedFileNameTURN = [resultsPathRoot species filesep 'BEDFILES' filesep 'TURN-' phenoName '.bed'];
bedFileIDTURN   = fopen(bedFileNameTURN,'w'); 

%MENT
bedFileNameMENT = [resultsPathRoot species filesep 'BEDFILES' filesep 'MENT-' phenoName '.bed'];
bedFileIDMENT   = fopen(bedFileNameMENT,'w'); 

%XENT
bedFileNameXENT = [resultsPathRoot species filesep 'BEDFILES' filesep 'XENT-' phenoName '.bed'];
bedFileIDXENT   = fopen(bedFileNameXENT,'w'); 

%PCC
%bedFileNamePCC = [resultsPathRoot species filesep 'BEDFILES' filesep 'PCC-' phenoName '.bed'];
%bedFileIDPCC   = fopen(bedFileNamePCC,'w'); 

% %ALPHA
% bedFileNameALPHA = [resultsPathRoot species filesep 'BEDFILES' filesep 'ALPHA-' phenoName '.bed'];
% bedFileIDALPHA   = fopen(bedFileNameALPHA,'w'); 
% 
% %BETA
% bedFileNameBETA = [resultsPathRoot species filesep 'BEDFILES' filesep 'BETA-' phenoName '.bed'];
% bedFileIDBETA   = fopen(bedFileNameBETA,'w'); 
% 
% %GAMMA
% bedFileNameGAMMA = [resultsPathRoot species filesep 'BEDFILES' filesep 'GAMMA-' phenoName '.bed'];
% bedFileIDGAMMA   = fopen(bedFileNameGAMMA,'w'); 

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write header info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% An
fprintf(bedFileIDAn,'%s%s\n',hardCodedOptionsAn,phenoName);

% Cn
fprintf(bedFileIDCn,'%s%s\n',hardCodedOptionsCn,phenoName);

% MNML
fprintf(bedFileIDMNML,'%s%s\n',hardCodedOptionsMNML,phenoName);

% Monostable
fprintf(bedFileIDMonostable,'%s%s\n',hardCodedOptionsMonostable,phenoName);

% Bistable
fprintf(bedFileIDBistable,'%s%s\n',hardCodedOptionsBistable,phenoName);

% NMEyClass
fprintf(bedFileIDNMEyClass,'%s%s\n',hardCodedOptionsNMEyClass,phenoName);

% NMEy
fprintf(bedFileIDNMEy,'%s%s\n',hardCodedOptionsNMEy,phenoName);

% MST
fprintf(bedFileIDMST,'%s%s\n',hardCodedOptionsMST,phenoName);

% MML
fprintf(bedFileIDMML,'%s%s\n',hardCodedOptionsMML,phenoName);

% ESI
fprintf(bedFileIDESI,'%s%s\n',hardCodedOptionsESI,phenoName);

% NCU
fprintf(bedFileIDNCU,'%s%s\n',hardCodedOptionsNCU,phenoName);

% NRDE
fprintf(bedFileIDNRDE,'%s%s\n',hardCodedOptionsNRDE,phenoName);

% TURN
fprintf(bedFileIDTURN,'%s%s\n',hardCodedOptionsTURN,phenoName);

% MENT
fprintf(bedFileIDMENT,'%s%s\n',hardCodedOptionsMENT,phenoName);

% XENT
fprintf(bedFileIDXENT,'%s%s\n',hardCodedOptionsXENT,phenoName);

% PCC
%fprintf(bedFileIDPCC,'%s%s\n',hardCodedOptionsPCC,phenoName);

% % ALPHA
% fprintf(bedFileIDALPHA,'%s%s\n',hardCodedOptionsALPHA,phenoName);
% 
% % BETA
% fprintf(bedFileIDBETA,'%s%s\n',hardCodedOptionsBETA,phenoName);
% 
% % GAMMA
% fprintf(bedFileIDGAMMA,'%s%s\n',hardCodedOptionsGAMMA,phenoName);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over chromosomes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


for chrNum=minChrNum:maxChrNum
    try
        chr_num_str = num2str(chrNum);
        chr_str     = ['chr' chr_num_str];

        %
        % Find last CpG site on chromosome (to determine how many regions to loop through)
        %
        CpGdata = [genomePathRoot 'genome' filesep species filesep 'CpGlocationChr' chr_num_str '.mat'];

        clear finalCpGloc CpGlocation %density           
        load(CpGdata,'finalCpGloc','CpGlocation');%,'density');

        %
        % find all start base-pairs for regions to be modeled on this chromsome
        %
        allStartBPs = int64(1):regionSize:int64(finalCpGloc); % only need to go to the last CpG site, not last BP

        %
        % Load relevant results files
        %
        clear mapObjData;
        load([resultsPathRoot species filesep chr_str filesep phenoName '_Analysis.mat'],'mapObjData');
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through all regions on chromosome
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        for regionNum = 1:length(allStartBPs)
            try
                %
                % find path of this region
                %
                startBP = allStartBPs(regionNum);
                endBP   = startBP+regionSize-int64(1);
                
                locationPathName = [chr_str '/bp' num2str(startBP) '-' num2str(endBP)];
                
                
                if isKey(mapObjData,locationPathName)
                    %
                    % load structs for region
                    %
                    regStruct   = mapObjData(locationPathName);
                    ProbsFromMC = regStruct.ProbsFromMC;  
                    isModeled   = regStruct.isModeled;
                    Nrq         = regStruct.Nrq;
                    NMEy        = regStruct.NMEy;
                    %NMEx        = regStruct.NMEx;
                    MedianVals  = double(regStruct.MedianVals);
                    MeanAn      = regStruct.MeanAn;
                    MeanCn      = regStruct.MeanCn;
                    MeanLevels  = regStruct.MeanLevels;
                    ESI         = regStruct.ESI;
                    NCU         = regStruct.NCU;
                    NRDE        = regStruct.NRDE;
                    TURN        = regStruct.TURN;
                    normEntReg  = double(regStruct.normEntReg);
                    MargEnt     = regStruct.MargEnt;
                    %PCC         = regStruct.NNcorr;
                    %thetabest   = regStruct.thetabest;
                    
%                     %
%                     % write out alpha beta and gamma
%                     %
%                     if ~isnan(thetabest(1))
%                         fprintf(bedFileIDALPHA,'%s %12u %12u %f\n',chr_str,int64(startBP-1),...
%                                                          int64(endBP-1),thetabest(1));
%                     end
%                     if ~isnan(thetabest(2))
%                         fprintf(bedFileIDBETA,'%s %12u %12u %f\n',chr_str,int64(startBP-1),...
%                                                          int64(endBP-1),thetabest(2));
%                     end
%                     if ~isnan(thetabest(3))
%                         fprintf(bedFileIDGAMMA,'%s %12u %12u %f\n',chr_str,int64(startBP-1),...
%                                                          int64(endBP-1),thetabest(3));
%                     end
                    
                    if ~isnan(normEntReg)
                        fprintf(bedFileIDXENT,'%s %12u %12u %f\n',chr_str,int64(startBP-1),...
                                                         int64(endBP-1),normEntReg);
                    end


                    %
                    % Compute regional quantities
                    %
                    
                    [lower_index,upper_index] = findSortedIndices(CpGlocation,startBP,endBP);
                    CpGlocInReg = CpGlocation(lower_index:upper_index);
                    numCpGinReg = length(CpGlocInReg);
                    
                    subRegCount = 0;
                    for offset = 0:subregionSize:(regionSize-1)
                        subRegCount = subRegCount + 1;
                        subStartBP  = startBP + offset;
                        subEndBP    = min(subStartBP + subregionSize - int64(1),finalCpGloc);
                     
                        numCpGsInSubReg = Nrq(subRegCount);
                        
                        if numCpGsInSubReg>0 && isModeled(subRegCount) % region modeled 
                            
                            [lower_index,upper_index] = findSortedIndices(CpGlocInReg,subStartBP,subEndBP);
                            if numCpGsInSubReg==1 && (lower_index==1||upper_index==numCpGinReg) % only 1 cpg site and its a boundary condition
                                continue; %move onto next subregion
                            end
                            
                            NMEyRegion    = NMEy(subRegCount);
                            yProbsRegion  = ProbsFromMC(subRegCount,:);
                            AnRegion      = MeanAn(subRegCount);
                            CnRegion      = MeanCn(subRegCount);
                            MNML          = MeanLevels(subRegCount);
                            ESIRegion     = ESI(subRegCount);
                            NCURegion     = NCU(subRegCount);
                            NRDEregion    = NRDE(subRegCount);
                            TURNregion    = TURN(subRegCount);
                            MargEntRegion = MargEnt(subRegCount);
                            %PCCregion    = PCC(subRegCount);
                            
                            %
                            % Print out An, Cn, MNML,ESI,NCU, NRDE
                            %
                            if AnRegion<inf
                                fprintf(bedFileIDAn,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),AnRegion);
                            end
                            if CnRegion<inf
                                fprintf(bedFileIDCn,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),CnRegion);
                            end
                            if MNML<inf
                                fprintf(bedFileIDMNML,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),MNML);
                            end
                            if  ESIRegion<inf
                                fprintf(bedFileIDESI,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),ESIRegion);
                            end
                            if NCURegion<inf
                                fprintf(bedFileIDNCU,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),NCURegion);
                            end
                            if NRDEregion<inf
                                fprintf(bedFileIDNRDE,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),NRDEregion);
                            end
                            if TURNregion<inf
                                fprintf(bedFileIDTURN,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),TURNregion);
                            end
                            if MargEntRegion<inf
                                fprintf(bedFileIDMENT,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),MargEntRegion);
                            end
%                             if PCCregion<inf
%                                 fprintf(bedFileIDPCC,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
%                                                        int64(subEndBP-1),PCCregion);
%                             end
                            
                            
                            
                            %
                            % Do NMEy classification
                            %
                            if NMEyRegion<inf
                                fprintf(bedFileIDNMEy,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),NMEyRegion);
                                if NMEyRegion >= 0.99 
                                    fprintf(bedFileIDNMEyClass,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),2);
                                elseif NMEyRegion >= 0.92 
                                    fprintf(bedFileIDNMEyClass,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),1);
                                elseif NMEyRegion > 0.444 
                                    fprintf(bedFileIDNMEyClass,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),0);
                                elseif NMEyRegion > 0.28
                                    fprintf(bedFileIDNMEyClass,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-1);
                                elseif NMEyRegion >=0    
                                    fprintf(bedFileIDNMEyClass,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-2);
                                end %end classification NMEy
                            end %end non-nan NMEy
                            
                            %
                            % Do methylation level classification
                            %
                            if  (numCpGsInSubReg > 1)
                                
                                unmethProb = yProbsRegion(1)+yProbsRegion(2);

                                fprintf(bedFileIDMST,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),MedianVals(subRegCount,1));
                                fprintf(bedFileIDMML,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),(MedianVals(subRegCount,1)/double(numCpGsInSubReg)));
                                
                                if unmethProb> (kappa/(1+kappa)) % monostable unmethylated

                                                                                                                           
                                    if yProbsRegion(1) >= kappa*(yProbsRegion(2)+yProbsRegion(3)+yProbsRegion(4)) % highly unmethylated
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),-2);
                                    else % unmethylated
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),-1);
                                    end
                                    
                                elseif unmethProb < (1/(1+kappa)) % monostable methylated
                                    
                                    
                                    if yProbsRegion(4) >= kappa*(yProbsRegion(1)+yProbsRegion(2)+yProbsRegion(3)) % highly methylated
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),2);
                                    else %  methylated
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),1);
                                    end
                                    
                                elseif ~(sum(isnan(yProbsRegion))>0) % mixture of phases (with no nan probabilities)
                                    
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),0); % print out 0 on METH track
        
                                    if (kappa*yProbsRegion(2) < yProbsRegion(1))&&(kappa*yProbsRegion(3) < yProbsRegion(4)) 
                                        %bistable
                                        
                                        fprintf(bedFileIDBistable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),3);
                                    elseif (kappa*yProbsRegion(1) < yProbsRegion(2))&&(kappa*yProbsRegion(4) < yProbsRegion(3))
                                        %mixed
                                                 
                                        fprintf(bedFileIDBistable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),1);
                                   
                                    elseif((yProbsRegion(1) < kappa*yProbsRegion(2))&&(yProbsRegion(2) < kappa*yProbsRegion(1))...
                                        &&(yProbsRegion(3) < kappa*yProbsRegion(4))&&(yProbsRegion(4) < kappa*yProbsRegion(3)))                                     
                                        %highly mixed
                                                  
                                        fprintf(bedFileIDBistable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),2);
                                         
                                    end % if well separated/mixed/other
                                    
                                end % if methylated/unmethylated/mixture
                                    
                                
                            else % numCpG==1
                                
                                unmethProb = yProbsRegion(1)+yProbsRegion(2);
                                
                                fprintf(bedFileIDMST,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),MedianVals(subRegCount,1));
                                fprintf(bedFileIDMML,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),MedianVals(subRegCount,1));
                                
                                if unmethProb> (kappa/(1+kappa)) % monostable unmethylated
                                              
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),-1);                                                                        
                                elseif unmethProb < (1/(1+kappa)) % monostable methylated
                                    
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),1);                                    
                                elseif ~(sum(isnan(yProbsRegion))>0) % mixed (with no nan probabilities)
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),0); 
                                    fprintf(bedFileIDBistable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                        int64(subEndBP-1),1);
                                end
                                
                            end % end min num of cpg for methylation classification                           
                                
                        end %region modeled
                        
                    end % loop over subregions
                end % end data present in region
            catch ME1
                disp(['Error at Chr' chr_num_str ':bp' num2str(startBP) '-' num2str(endBP) ' with stack:'])
                ME1.stack.file
                ME1.stack.name
                ME1.stack.line
                ME1.identifier
            end
        end %loop over regions
    catch ME2
        disp(['Error at Chr ' num2str(chrNum) ' with stack:'])
        ME2.stack.file
        ME2.stack.name
        ME2.stack.line
        ME2.identifier
    end
end % loop over chromosomes

%
% Close all files currently open for writing
%

fclose(bedFileIDAn);
fclose(bedFileIDCn);
fclose(bedFileIDMNML);
fclose(bedFileIDMonostable);
fclose(bedFileIDBistable);
fclose(bedFileIDNMEyClass);
fclose(bedFileIDNMEy);
fclose(bedFileIDMST);
fclose(bedFileIDMML);
fclose(bedFileIDESI);
fclose(bedFileIDNCU);
fclose(bedFileIDNRDE);
fclose(bedFileIDTURN);
fclose(bedFileIDMENT);
fclose(bedFileIDXENT);
%fclose(bedFileIDPCC);
% fclose(bedFileIDALPHA);
% fclose(bedFileIDBETA);
% fclose(bedFileIDGAMMA);

end