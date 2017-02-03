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
%%%%%%%%                    processMatrix.m                        %%%%%%%%
%%%%%%%%          Code written by: W. Garrett Jenkinson            %%%%%%%%
%%%%%%%%               Last Modified: 12/01/2016                   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a data matrix and breaks non-contiguous 
% observations of CpG sites into independent reads. It also provides 
% the start and end indices of each read, allowing the data matrix 
% to be parsed more quickly in downstream applications. 
%
% USAGE:
%
% [newMatrix,CpGstart,CpGend] = processMatrix(dataMat)
%
% INPUT:
%
% dataMat
%           A matrix with each row corresponding to a methylation read 
%           and each column corresponding to a CpG site. The elements 
%           of this matrix take values -1, 0, or 1 which indicate that 
%           a given read does not observe a CpG site (-1), it observes 
%           lack of methylation at a CpG site (0), or it observes 
%           methylation at a CpG site (1). 
%
% OUTPUTS:
%
% newMatrix
%           A matrix with each row having only contiguous observations of
%           CpG sites. The resulting matrix has the same number of
%           observations as dataMat, but with the new requirement that
%           reads have only contiguous CpG sites observed (and thus a 
%           read in dataMat that has a gap in the observations will be 
%           broken into multiple reads in newMatrix).
%
% CpGstart
%           A vector with as many elements as the number of rows in 
%           newMatrix. Each element of this vector contains the index 
%           of the first observation in the corresponding row of 
%           newMatrix.
%
% CpGend
%           A vector with as many elements as the number of rows in 
%           newMatrix. Each element of this vector contains the index 
%           of the last observation in the corresponding row of 
%           newMatrix.
%

function  [newMatrix,CpGstart,CpGend] = processMatrix(dataMat)

% Remove empty rows of data. 

numEmptyRows = sum(sum(dataMat>-1,2)==0);
if numEmptyRows>0
    dataMat = dataMat(sum(dataMat>-1,2)>0,:);
    disp(['WARNING: Data matrix contains ' num2str(numEmptyRows) ' empty rows'])
end

% Get number of reads in data matrix and initialize output variables.

origNumReads = size(dataMat,1); % Original number of reads.

% Preallocate space. If more space is needed it is handled in a try-catch
% block in the loop below.

newMatrix = -1*ones(2*origNumReads,size(dataMat,2)); 
CpGstart  = zeros(2*origNumReads,1);
CpGend    = zeros(2*origNumReads,1);

% Preprocess read data to deal with gaps in the middle of a read.

readCountNew = 0; % Count for the number of reads in newMatrix.

for readNum = 1:origNumReads % Process each read in dataMat individually.
    
    % Location of previously observed CpG site. Set to -2 since 
    % we have no valid previous observations.
    prevObservation = -2; 
    
    % Flag to determine if we are at the first observation on the 
    % current dataMat read.
    firstObservationOnRead = 1;
    
    for observation = find(dataMat(readNum,:)>-1) 
              % Loop through each relevant observation on this dataMat
              % read.
        if (observation-prevObservation)==1 % Still in same contiguous portion.

            newMatrix(readCountNew,observation) = dataMat(readNum,observation);
              % Take observation from dataMat and put it in newMatrix.
            prevObservation = observation; % Set up for next iteration of loop.

        else % New contiguous portion.

            readCountNew = readCountNew+1; 
                 % New contiguous portion means another count for
                 % newMatrix.
            try
                % Store observation in newMatrix and record start of 
                % new CpG observation.
                newMatrix(readCountNew,observation) = dataMat(readNum,observation);
                CpGstart(readCountNew)              = observation; 

            catch % Exceed newMatrix dimension, increase matrix size first 
                  % and try again.

                newMatrix = [newMatrix;-1*ones(2*origNumReads,...
                                             size(dataMat,2))]; %#ok<AGROW>
                CpGstart  = [CpGstart;zeros(2*origNumReads,1)]; %#ok<AGROW>
                CpGend    = [CpGend;zeros(2*origNumReads,1)];   %#ok<AGROW>
                
                newMatrix(readCountNew,observation) = dataMat(readNum,observation);
                CpGstart(readCountNew)              = observation; 

            end % End try catch.
            
            if ~firstObservationOnRead
                CpGend(readCountNew-1) = prevObservation; 
                                          % Record the end of previous
                                          % read.
            else
                firstObservationOnRead=0; % Next iteration of loop will 
                                          % no longer be the first 
                                          % observation on the read.
            end

            prevObservation = observation; % Set up for next iteration of loop.
         
        end % End test for contiguous portion.
        
    end % End loop through observed CpGs in current read.
    
    CpGend(readCountNew)=observation; % The index of the last observation
                                      % will be the end of the current
                                      % read.
        
end

% Resize final output to only have relevant information (i.e., deal with
% the unused preallocated space).

newMatrix = newMatrix(1:readCountNew,:);
CpGstart  = CpGstart(1:readCountNew);
CpGend    = CpGend(1:readCountNew);