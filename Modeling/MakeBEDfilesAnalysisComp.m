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
%%%%%%%%%%%             Last Modified: 05/09/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function makes the relevant bed files to compare two phenotypes. It
% assumes that statistical estimation has been completed using
% EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
% with AnalysisForChr.m/MergeAnalysis.m for both phenotypes.
%
% Example Default Usage:
% MakeBEDfilesAnalysisComp(phenoNames)
%
% Example Usage Modifying optional parameter "maxChrNum" as name-value pair:
% MakeBEDfilesAnalysisComp(phenoNames,'maxChrNum',5)
%
% This file takes as mandatroy inputs (which must appear first and in the 
% following order):
%
% phenoNames
%               A 2x1 cell with strings of the name of the 
%               phenotypes that have already been computed by
%               EstParamsForChr.m
% 
% This file takes as optional inputs, specified in any order after the 
% mandatory inputs as name-value pairs; i.e.,
% MakeBEDfilesAnalysisComp(...,'inputName',inputValue):
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

function MakeBEDfilesAnalysisComp(phenoNames,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'phenoNames')
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
addParameter(p,'minNumCpG',2,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','nonnegative','scalar'}))
addParameter(p,'zProbThresh',0.55,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','nonnegative','scalar'}))
addParameter(p,'boundariesZ',[-1,-0.55,-0.1,0.1,0.55,1],...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','nonnegative','scalar','size',[1,6]}))
addParameter(p,'boundariesDER',[-1,-.5,-.3,-.05,.05,.3,.5,1],...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','nonnegative','scalar','size',[1,8]}))

           
parse(p,phenoNames,varargin{:})

species         = p.Results.species;
minChrNum       = p.Results.minChrNum;
maxChrNum       = p.Results.maxChrNum;
resultsPathRoot = p.Results.resultsPathRoot;
genomePathRoot  = p.Results.genomePathRoot;
regionSize      = p.Results.regionSize;
subregionSize   = p.Results.subregionSize;
minNumCpG       = p.Results.minNumCpG;
zProbThresh     = p.Results.zProbThresh;
boundariesZ     = p.Results.boundariesZ;
boundariesDER   = p.Results.boundariesDER;


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

%Monostable
hardCodedOptionsMonostable = 'track type=bedGraph visibility=full autoScale=off viewLimits=-3.0:3.0 name=DMU-';

%NMEy
hardCodedOptionsNMEy = 'track type=bedGraph visibility=full autoScale=off viewLimits=-1.0:1.0 name=dNME-';

%dENTR
hardCodedOptionsENTR = 'track type=bedGraph visibility=full autoScale=off viewLimits=-3.0:3.0 name=DEU-';

%MST
hardCodedOptionsMST = 'track type=bedGraph visibility=full autoScale=on name=dMDS-';

%MML
hardCodedOptionsMML = 'track type=bedGraph visibility=full autoScale=off viewLimits=-1.0:1.0 name=dMDL-';

%MNML
hardCodedOptionsMNML = 'track type=bedGraph visibility=full autoScale=off viewLimits=-1.0:1.0 name=dMNL-';

%ESI
hardCodedOptionsESI = 'track type=bedGraph visibility=full autoScale=on name=dESI-';

%dNCU
hardCodedOptionsNCU = 'track type=bedGraph visibility=full autoScale=off viewLimits=-1.0:1.0 name=dCAP-';

%dNRDE
hardCodedOptionsNRDE = 'track type=bedGraph visibility=full autoScale=on name=dRDE-';

%Dpq
hardCodedOptionsJS = 'track type=bedGraph visibility=full autoScale=off viewLimits=0.0:1.0 name=JSD-';


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

%Monostable
bedFileNameMonostable = [resultsPathRoot species filesep 'BEDFILES' filesep 'DMU-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDMonostable   = fopen(bedFileNameMonostable,'w'); 

%NMEy
bedFileNameNMEy = [resultsPathRoot species filesep 'BEDFILES' filesep 'dNME-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDNMEy   = fopen(bedFileNameNMEy,'w'); 

%ENTR
bedFileNameENTR = [resultsPathRoot species filesep 'BEDFILES' filesep 'DEU-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDENTR   = fopen(bedFileNameENTR,'w'); 

%MST
bedFileNameMST = [resultsPathRoot species filesep 'BEDFILES' filesep 'dMDS-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDMST   = fopen(bedFileNameMST,'w'); 

%MML
bedFileNameMML = [resultsPathRoot species filesep 'BEDFILES' filesep 'dMDL-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDMML   = fopen(bedFileNameMML,'w'); 

%MNML
bedFileNameMNML = [resultsPathRoot species filesep 'BEDFILES' filesep 'dMNL-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDMNML   = fopen(bedFileNameMNML,'w'); 

%ESI
bedFileNameESI = [resultsPathRoot species filesep 'BEDFILES' filesep 'dESI-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDESI   = fopen(bedFileNameESI,'w'); 

%NCU
bedFileNameNCU = [resultsPathRoot species filesep 'BEDFILES' filesep 'dCAP-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDNCU   = fopen(bedFileNameNCU,'w'); 

%NRDE
bedFileNameNRDE = [resultsPathRoot species filesep 'BEDFILES' filesep 'dRDE-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDNRDE   = fopen(bedFileNameNRDE,'w'); 

%JS dist
bedFileNameJS = [resultsPathRoot species filesep 'BEDFILES' filesep 'JSD-' phenoNames{1} '-VS-' phenoNames{2} '.bed'];
bedFileIDJS   = fopen(bedFileNameJS,'w');

    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write header info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Monostable
fprintf(bedFileIDMonostable,'%s%s-VS-%s\n',hardCodedOptionsMonostable,phenoNames{1},phenoNames{2});

% NMEy
fprintf(bedFileIDNMEy,'%s%s-VS-%s\n',hardCodedOptionsNMEy,phenoNames{1},phenoNames{2});

% ENTR
fprintf(bedFileIDENTR,'%s%s-VS-%s\n',hardCodedOptionsENTR,phenoNames{1},phenoNames{2});

% MST
fprintf(bedFileIDMST,'%s%s-VS-%s\n',hardCodedOptionsMST,phenoNames{1},phenoNames{2});

% MML
fprintf(bedFileIDMML,'%s%s-VS-%s\n',hardCodedOptionsMML,phenoNames{1},phenoNames{2});

% MNML
fprintf(bedFileIDMNML,'%s%s-VS-%s\n',hardCodedOptionsMNML,phenoNames{1},phenoNames{2});

% ESI
fprintf(bedFileIDESI,'%s%s-VS-%s\n',hardCodedOptionsESI,phenoNames{1},phenoNames{2});

% NCU
fprintf(bedFileIDNCU,'%s%s-VS-%s\n',hardCodedOptionsNCU,phenoNames{1},phenoNames{2});

% NRDE
fprintf(bedFileIDNRDE,'%s%s-VS-%s\n',hardCodedOptionsNRDE,phenoNames{1},phenoNames{2});

% JS dist
fprintf(bedFileIDJS,'%s%s-VS-%s\n',hardCodedOptionsJS,phenoNames{1},phenoNames{2});


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

        clear finalCpGloc CpGlocation          
        load(CpGdata,'finalCpGloc','CpGlocation');%,'density');

        %
        % find all start base-pairs for regions to be modeled on this chromsome
        %
        allStartBPs = int64(1):regionSize:int64(finalCpGloc); % only need to go to the last CpG site, not last BP

        %
        % Load relevant results files
        %
        clear mapObjData;
        load([resultsPathRoot species filesep chr_str filesep phenoNames{1} '_Analysis.mat'],'mapObjData');
        mapObjData1 = mapObjData;
        clear mapObjData;
        load([resultsPathRoot species filesep chr_str filesep phenoNames{2} '_Analysis.mat'],'mapObjData');
        mapObjData2 = mapObjData;
        clear mapObjData;
        
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
                
                
                if isKey(mapObjData1,locationPathName)&&isKey(mapObjData2,locationPathName)
                    %
                    % load structs for region
                    %
                    regStruct1   = mapObjData1(locationPathName);
                    isModeled1   = regStruct1.isModeled;
                    NMEy1        = regStruct1.NMEy;
                    fullProbs1   = regStruct1.fullProbs;
                    %MedianVals1  = double(regStruct1.MedianVals);
                    Nrq1         = regStruct1.Nrq;
                    ESI1         = regStruct1.ESI;
                    NCU1         = regStruct1.NCU;
                    NRDE1        = regStruct1.NRDE;
                    
                    regStruct2   = mapObjData2(locationPathName); 
                    isModeled2   = regStruct2.isModeled;
                    NMEy2        = regStruct2.NMEy;
                    fullProbs2   = regStruct2.fullProbs;
                    %MedianVals2  = double(regStruct2.MedianVals);
                    Nrq2         = regStruct2.Nrq;
                    ESI2         = regStruct2.ESI;
                    NCU2         = regStruct2.NCU;
                    NRDE2        = regStruct2.NRDE;
                    
                    if sum(Nrq1~=Nrq2)
                        disp('ERROR: N1 does not equal N2');
                        ME = MException('VerifyOutput:OutOfBounds', ...
                            'Results are outside the allowable limits');
                        throw(ME);
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
                        subStartBP = startBP+offset;
                        subEndBP   = min(subStartBP + subregionSize - int64(1),finalCpGloc);
                        
                        numCpGsInSubReg = Nrq1(subRegCount);
                        
                        if numCpGsInSubReg>0 && (isModeled1(subRegCount)) && (isModeled2(subRegCount)) % region modeled 
                            
                            [lower_index,upper_index] = findSortedIndices(CpGlocInReg,subStartBP,subEndBP);
                            if numCpGsInSubReg==1 && (lower_index==1||upper_index==numCpGinReg) % only 1 cpg site and its a boundary condition
                                continue; %move onto next subregion
                            end
                            
                            
                            NMEyRegion1   = NMEy1(subRegCount);
                            NMEyRegion2   = NMEy2(subRegCount);
                            dNMEyRegion   = NMEyRegion1-NMEyRegion2;
                            if dNMEyRegion<inf %no nans or infs
                                fprintf(bedFileIDNMEy,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),dNMEyRegion);
                            end
                            
                                                                           
                            dNCU  = NCU1(subRegCount)-NCU2(subRegCount);
                            if dNCU<inf %no nans or infs
                                fprintf(bedFileIDNCU,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),dNCU);
                            end
                                        
                            dNRDE  = NRDE1(subRegCount)-NRDE2(subRegCount);
                            if dNRDE<inf %no nans or infs
                                fprintf(bedFileIDNRDE,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),dNRDE);
                            end
                            
                            dESI  = ESI1(subRegCount)-ESI2(subRegCount);
                            if dESI<inf %strongly hyposensitive
                                fprintf(bedFileIDESI,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),dESI);
                            end
                            
                                               
                            %
                            % Do NMEy classification with:
                            %  -1 <= dNME <= -0.5  %high
                            % -0.5 <  dNME <= -0.3 %moderate
                            % -0.3 <  dNME <= -0.05 % weak
                            % -0.05 <  dNME <   0.05 %iso
                            %  0.05 <= dNME <   0.3  % weak
                            %  0.3 <= dNME <   0.5  %moderate
                            %  0.5 <= dNME <=  1   %high
                            %
                            
                            
                            if dNMEyRegion >= boundariesDER(7) %high
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),3);
                            elseif dNMEyRegion >= boundariesDER(6) %moderate
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),2);
                            elseif dNMEyRegion >= boundariesDER(5) %weak
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),1);
                            elseif dNMEyRegion > boundariesDER(4) %iso
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),0);
                            elseif dNMEyRegion > boundariesDER(3) %weak
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),-1);
                            elseif dNMEyRegion > boundariesDER(2)   %moderate
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),-2);
                            elseif dNMEyRegion >=boundariesDER(1)   %high
                                fprintf(bedFileIDENTR,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),-3);
                            end %end classification NMEy
                            
                            %
                            % Compute distribution of Y1-Y2 and its mean
                            %
                            pY1 = full(fullProbs1(subRegCount,1:(numCpGsInSubReg+1)));
                            pY2 = full(fullProbs2(subRegCount,1:(numCpGsInSubReg+1)));
                            
                            pY1minusY2 = conv(pY1,pY2(end:-1:1));
                            zVals      = -1:(1/double(numCpGsInSubReg)):1;
                            
                            %Find mean
                            dMNML = dot(pY1minusY2,zVals);
                            if dMNML<inf %no nans or infinities
                                fprintf(bedFileIDMNML,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),dMNML);
                            end
                                               
                            % Find median
                            med_index = find(cumsum(pY1minusY2)>=0.5,1,'first');
                            dMDML     = zVals(med_index);
                            dMDMS     = dMDML*numCpGsInSubReg;
                            if dMDMS<inf %no nans or infinities
                                fprintf(bedFileIDMST,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),dMDMS);
                                fprintf(bedFileIDMML,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                   int64(subEndBP-1),dMDML);                                              
                            end    
                            
                            %
                            % Compute JS distance
                            %
                            
                            pY1=pY1(:);%ensures column vector
                            pY2=pY2(:);%ensures column vector
                           
                            % JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
                            % where M=(P+M)/2
                            % and KL(P|M)=sum_i P_i log2(P_i/M_i)
                            Dpq = sqrt( sum( pY1(pY1>0).*log2(2.*pY1(pY1>0)./(pY1(pY1>0)+pY2(pY1>0))) )/2 ...
                                       +sum( pY2(pY2>0).*log2(2.*pY2(pY2>0)./(pY1(pY2>0)+pY2(pY2>0))) )/2 );
                            
                            if Dpq<inf % ensure no nans or infinities         
                                fprintf(bedFileIDJS,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),Dpq);
                            end
                                                                           
                            %
                            % Compute coarse probabilities q
                            %                              
                            
                            zProbsRegion    = zeros(5,1);
                            zProbsRegion(1) = sum(pY1minusY2((boundariesZ(1)<=zVals)&(zVals<=boundariesZ(2)))); %   -1 <= Z <= -0.8
                            zProbsRegion(2) = sum(pY1minusY2((boundariesZ(2)<zVals)&(zVals<=boundariesZ(3))));  % -0.8 <  Z <= -0.1
                            zProbsRegion(3) = sum(pY1minusY2((boundariesZ(3)<zVals)&(zVals<boundariesZ(4))));   % -0.1 <  Z <   0.1
                            zProbsRegion(4) = sum(pY1minusY2((boundariesZ(4)<=zVals)&(zVals<boundariesZ(5))));  %  0.1 <= Z <   0.8
                            zProbsRegion(5) = sum(pY1minusY2((boundariesZ(5)<=zVals)&(zVals<=boundariesZ(6)))); %  0.8 <= Z <=  1

                            if numCpGsInSubReg>=minNumCpG % enough to do normal classification
                                if zProbsRegion(3) > zProbThresh
                                    %ISO-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),0);
                                elseif (zProbsRegion(4)+zProbsRegion(5)) > zProbThresh
                                    if zProbsRegion(5) > zProbThresh
                                        %Ph-1 STRONGLY HYPER-METHYLATED
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),3);
                                    else
                                        %Ph-1 HYPER-METHYLATED
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),2);
                                    end

                                elseif (zProbsRegion(1)+zProbsRegion(2)) > zProbThresh
                                    if zProbsRegion(1) > zProbThresh
                                        %Ph-1 STRONGLY HYPO-METHYLATED
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-3);
                                    else
                                        %Ph-1 HYPO-METHYLATED
                                        fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-2);
                                    end

                                elseif (zProbsRegion(4)+zProbsRegion(5))/(1-zProbsRegion(3)) > zProbThresh% weakly hypometh
                                    %Ph-1 WEAKLY HYPER-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),1);
                                elseif (zProbsRegion(1)+zProbsRegion(2))/(1-zProbsRegion(3)) > zProbThresh % weakly hypermeth
                                    %Ph-1 WEAKLY HYPO-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-1);
                                %else %No classification
                                  %N/C

                                end
                            else % no "highly" classification allowed
                                if zProbsRegion(3) > zProbThresh
                                    %ISO-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),0);
                                elseif (zProbsRegion(4)+zProbsRegion(5)) > zProbThresh
                                    %Ph-1 HYPER-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),2);
                                elseif (zProbsRegion(1)+zProbsRegion(2)) > zProbThresh
                                     %Ph-1 HYPO-METHYLATED
                                     fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-2);
                                elseif (zProbsRegion(4)+zProbsRegion(5))/(1-zProbsRegion(3)) > zProbThresh% weakly hypometh
                                    %Ph-1 WEAKLY HYPER-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),1);
                                 elseif (zProbsRegion(1)+zProbsRegion(2))/(1-zProbsRegion(3)) > zProbThresh % weakly hypermeth
                                    %Ph-1 WEAKLY HYPO-METHYLATED
                                    fprintf(bedFileIDMonostable,'%s %12u %12u %f\n',chr_str,int64(subStartBP-1),...
                                                       int64(subEndBP-1),-1);
                                end
                                    
                            end % if numCpGsInSubReg>=minNumCpG, else
                                                           
                        end %region modeled
                        
                        
                    end % loop over subregions
                end % data in both regions
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
end %loop over chromosomes

%
% Close all files currently open for writing
%

fclose(bedFileIDMonostable);
fclose(bedFileIDNMEy);
fclose(bedFileIDENTR);
fclose(bedFileIDMST);
fclose(bedFileIDMML);
fclose(bedFileIDMNML);
fclose(bedFileIDESI);
fclose(bedFileIDNCU);
fclose(bedFileIDNRDE);
fclose(bedFileIDJS);

end