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
%%%%%%%%%%%             Last Modified: 09/24/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function does PostProcessing of bed files to find MBs. 
% It assumes that statistical estimation has been completed using
% EstParamsForChr.m/MergeEstParams.m and that these results were analyzed 
% with AnalysisForChr.m/MergeAnalysis.m  and bed files generated with
% MakeBEDfilesAnalysis.m.
%
% Example usage: 
% MethBlocks( 'colonnormal', 'Human')
%
% Inputs:
%
% phenoName
%           The name of the phenotype for which bed files have been made.
%
% species
%           The species of the phenotype. For now only 'Human' and 'Mouse'
%           are accepted, since the algorithms require specification
%           of chromosome lengths.
%
function MethBlocks( phenoName, species ,varargin)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse values passed as inputs to the fuction and validate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = inputParser;

addRequired(p,'phenoName')
addRequired(p,'species')
addParameter(p,'resultsPathRoot',['.' filesep 'results' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'genomePathRoot',['..' filesep 'ParseBAMfile' filesep],...
               @(x)validateattributes(x,{'char'},{'nonempty'}))
addParameter(p,'subRegionSize',150,...
               @(x)validateattributes(x,{'numeric'},{'nonempty','integer','positive','scalar'}))

parse(p,phenoName, species,varargin{:})

resultsPathRoot = p.Results.resultsPathRoot;
genomePathRoot  = p.Results.genomePathRoot;
subRegionSize   = p.Results.subRegionSize; 


%
% Manual checks/corrections of inputs
%

if genomePathRoot(end)~=filesep
    genomePathRoot=[genomePathRoot filesep];
end
if resultsPathRoot(end)~=filesep
    resultsPathRoot=[resultsPathRoot filesep];
end

trackPath = [resultsPathRoot species '/BEDFILES/'];



%
% Initialization
%

threshVal      = 3/4; % percentage of disordered regions required to make a disorder region 
numRegInWindow = 500;

threshNumRegModeled = 0.075; 



%
% Load data from files
%
bedFileName = [trackPath 'METH-' phenoName '.bed'];

FID = fopen(bedFileName);
formatSpec='chr%u %u %u %f';
bedData = textscan(FID,formatSpec,'HeaderLines',1);

fclose(FID);



%
% Open file for writing results
%

bedFileIDfiltered = fopen([trackPath 'MB-'  phenoName  '.bed'],'w'); %  ...
              %'_window' num2str(numRegInWindow) '_thresh' num2str(threshVal) '.bed'],'w');

hardCodedOptions = ['track  type=bedGraph visibility=full maxHeightPixels=100:15:11 autoScale=off viewLimits=-1.0:1.0 color=190,0,0 altColor=0,140,0 name=MB-' phenoName ];%...
            %'_window' num2str(numRegInWindow) '_thresh' num2str(threshVal) ];
fprintf(bedFileIDfiltered,'%s\n',hardCodedOptions);


%
% Loop through regions in bed file
%

chrList = unique(bedData{1});

for bedChr = chrList'
    % get data on the current chromosome
    bedStartBPcurrChr = bedData{2}(bedData{1}==bedChr);
    bedValsCurrChr    = bedData{4}(bedData{1}==bedChr);
    disorderedCurrChr = bedValsCurrChr>0;
    orderedCurrChr    = bedValsCurrChr<0;
    
    minWindStartBPchr = max(1,bedStartBPcurrChr(1)-(numRegInWindow-1)*subRegionSize);
    maxWindStartBPchr = bedStartBPcurrChr(end)+(numRegInWindow-1)*subRegionSize;
    
    % initialize 
    prevInsideDisRegion  = 0;
    prevDisRegionStartBP = inf;
    prevDisRegionEndBP   = 0;
    prevInsideOrdRegion  = 0;
    prevOrdRegionStartBP = inf;
    prevOrdRegionEndBP   = 0;
    mostRecentConflict   = 0;
    currRegionStartBP    = inf;
    currRegionEndBP      = 0;
    
    % load last CpG location on chromosome
    CpGdata = [genomePathRoot 'genome' filesep species filesep 'CpGlocationChr' num2str(bedChr) '.mat'];
    load(CpGdata,'chrLength');

    for currMinWindowStartBP = minWindStartBPchr:subRegionSize:maxWindStartBPchr
        currMaxWindowStartBP = currMinWindowStartBP +(numRegInWindow-1)*subRegionSize;
        
        % compute the number of disordered units within window
        [lower_index,upper_index] = findSortedIndices(bedStartBPcurrChr,currMinWindowStartBP,currMaxWindowStartBP);
        disorderCounts = sum(disorderedCurrChr(lower_index:upper_index));
        orderCounts    = sum(orderedCurrChr(lower_index:upper_index));
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check it is a disordered window or ordered window
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        currInsideDisRegion=0;
        currInsideOrdRegion=0;
        if (orderCounts+disorderCounts)>=threshNumRegModeled*numRegInWindow % this has enough modeled
            if disorderCounts>threshVal*(orderCounts+disorderCounts)% This is a disordered region
                % Update Region information 
                currInsideDisRegion  = 1;
                currRegionStartBP = currMinWindowStartBP;
                currRegionEndBP   = currMaxWindowStartBP+subRegionSize-1;
            elseif orderCounts>threshVal*(orderCounts+disorderCounts)% this is a ordered region
                % Update Region information 
                currInsideOrdRegion  = 1;
                currRegionStartBP = currMinWindowStartBP;
                currRegionEndBP   = currMaxWindowStartBP+subRegionSize-1;
            end
        end
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure out what to do based on current and previous region status
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if (currRegionStartBP>mostRecentConflict) &&... % current region is not a conflict region 
            (currInsideDisRegion||currInsideOrdRegion)
            if prevInsideDisRegion
                if currInsideDisRegion % had disorder previously, and now again
                    if prevDisRegionEndBP<currRegionStartBP 
                        % current region does not overlap previous region
                        % writeout previous region and update for next iteration
                        fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevDisRegionStartBP-1),int64(prevDisRegionEndBP-1),+1);

                        prevDisRegionStartBP = currRegionStartBP;
                        prevDisRegionEndBP   = currRegionEndBP;
                    else
                        % current region overlaps previous region, extend end of
                        % previous region
                        prevDisRegionEndBP = currRegionEndBP;
                    end
                elseif currInsideOrdRegion
                    % Disordered region ends and ordered region begins or
                    % conflict region
                    prevInsideDisRegion = 0;
                    prevInsideOrdRegion = 1;
                    mostRecentConflict  = prevDisRegionEndBP;
                    
                    %writeout previous disordered region, not including overlap
                    %with the current ordered region
                    fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevDisRegionStartBP-1),...
                             min(int64(prevDisRegionEndBP-1),int64(currRegionStartBP-1)),+1);

                    % setup current ordered region, without including overlap
                    % of previous region
                    prevOrdRegionStartBP = max(prevDisRegionEndBP,currRegionStartBP);
                    prevOrdRegionEndBP   = currRegionEndBP;

                    % reset previous disoredered region 
                    prevDisRegionStartBP = inf;
                    prevDisRegionEndBP   = 0; 
                %else % not currently inside ordered or disordered region
                    % do nothing so we comment this out
                end
            elseif prevInsideOrdRegion
                if currInsideDisRegion
                    % Ordered region ends and disordered region begins or
                    % conflict region
                    prevInsideDisRegion = 1;
                    prevInsideOrdRegion = 0;
                    mostRecentConflict  = prevOrdRegionEndBP;

                    %writeout previous disordered region, not including overlap
                    %with the current ordered region
                    fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevOrdRegionStartBP-1),...
                             min(int64(prevOrdRegionEndBP-1),int64(currRegionStartBP-1)),-1);


                    % setup current disordered region, without including overlap
                    % of previous region
                    prevDisRegionStartBP = max(prevOrdRegionEndBP,currRegionStartBP);
                    prevDisRegionEndBP   = currRegionEndBP;

                    % reset previous oredered region 
                    prevOrdRegionStartBP = inf;
                    prevOrdRegionEndBP   = 0; 

                elseif currInsideOrdRegion
                    % had order before, and now again
                    if prevOrdRegionEndBP<currRegionStartBP 
                        % current region does not overlap previous region
                        % writeout previous region and update for next iteration
                        fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevOrdRegionStartBP-1),int64(prevOrdRegionEndBP-1),-1);

                        prevOrdRegionStartBP = currRegionStartBP;
                        prevOrdRegionEndBP   = currRegionEndBP;
                    else
                        % current region overlaps previous region, extend end of
                        % previous region
                        prevOrdRegionEndBP = currRegionEndBP;
                    end
                %else % not currently inside ordered or disordered region
                    %  do nothing, so we comment this case out
                end
            else % previously not in ordered or disordered region
                if currInsideDisRegion
                    % new disordered region begins
                    prevInsideDisRegion  = 1;
                    prevDisRegionStartBP = currRegionStartBP;
                    prevDisRegionEndBP   = currRegionEndBP;
                elseif currInsideOrdRegion
                    % new ordered region begins
                    prevInsideOrdRegion  = 1;
                    prevOrdRegionStartBP = currRegionStartBP;
                    prevOrdRegionEndBP   = currRegionEndBP;
                %else %: not currently inside ordered or disordered region
                    % Nothing happens so we ignore this case
                end
            end % end case handling for previous/current values
        end % current region not a conflict region
        
    end % loop over start BPs
    
    %
    % If we reached the end of chromosome, make sure we writeout the last
    % ordered/disordered region we had (if not already done)
    %
    if prevInsideDisRegion
        if int64(prevDisRegionStartBP-1)<min(chrLength,int64(prevDisRegionEndBP-1))
            fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevDisRegionStartBP-1),...
                                 min(chrLength,int64(prevDisRegionEndBP-1)),+1);
        end
    elseif prevInsideOrdRegion 
        if int64(prevOrdRegionStartBP-1)<min(chrLength,int64(prevOrdRegionEndBP-1))
            fprintf(bedFileIDfiltered,'%s %12u %12u %f\n',['chr' num2str(bedChr)],int64(prevOrdRegionStartBP-1),...
                                 min(chrLength,int64(prevOrdRegionEndBP-1)),-1);
        end
    end
    
end

fclose(bedFileIDfiltered);

end
