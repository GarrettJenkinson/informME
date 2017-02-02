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
%%%%%%%%%%%             Last Modified: 05/15/2015              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function takes a data matrix and breaks apart non-contiguous
% observations of CpG sites into independent reads. It also provides the
% start and end index of each read, allowing the data matrix to be parsed
% more quickly in downstream applications. 
%
% [newMatrix, CpGstart, CpGend] = processMatrix( dataMat );
%
% Takes as input:
%
% dataMat
%           A matrix with each row corresponding to a read and each column
%           corresponding to a CpG site. It takes values of -1, 0, 1
%           which indicate that a given read does not observe a CpG site,
%           observes a lack of methylation at a CpG site, or observes
%           methylation at a CpG site, respectively.
%
% Produces as outputs:
%
% newMatrix
%           A matrix with each row having only contiguous observations of
%           CpG sites. The resulting matrix has the same number of
%           observations as dataMat, but with the new requirement that
%           reads only have contiguous CpG sites observed (and thus a read
%           in dataMat that had a gap in the observations will be broken up
%           into multiple reads in newMatrix).
%
% CpGstart
%           A vector with as many elements as rows of newMatrix. Each
%           element contains the index of the first observation in the
%           corresponding row of newMatrix.
%
% CpGend
%           A vector with as many elements as rows of newMatrix. Each
%           element contains the index of the last observation in the
%           corresponding row of newMatrix.

function  [newMatrix, CpGstart, CpGend] = processMatrix( dataMat )



%
% First remove any empty rows of data (should not be any, defensive programming)
%
numEmptyRows = sum(sum(dataMat>-1,2)==0);
if numEmptyRows>0
    dataMat = dataMat(sum(dataMat>-1,2)>0,:);
    display(['WARNING: Data matrix contains ' num2str(numEmptyRows) ' empty rows'])
end

%
% Now get number of reads in data matrix and initialize output vars
%
origNumReads = size(dataMat,1); % original number of reads


% preallocate space...if more space is needed it is handled in a try-catch
% block in the loop below
newMatrix = -1*ones(2*origNumReads,size(dataMat,2)); 
CpGstart  = zeros(2*origNumReads,1);
CpGend    = zeros(2*origNumReads,1);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess read data to deal with gaps in the middle of a read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

readCountNew = 0; % the count for the number of reads in newMatrix

for readNum = 1:origNumReads %process each of the reads in dataMat individually
    
    % location of previously observed CpG site...set to -2 since we have no valid previous observations
    prevObservation = -2; 
    
    % flag to determine if we are at the first observation on the current dataMat read
    firstObservationOnRead = 1;
    
    for observation = find(dataMat(readNum,:)>-1) % loop through each relevant observation on this dataMat read
        if (observation-prevObservation)==1 % still in same contiguous portion
            
            % take observation from dataMat and put it in newMatrix
            newMatrix(readCountNew,observation) = dataMat(readNum,observation); 
            
            prevObservation = observation; % set up for next iteration of loop
            
        else % new contiguous portion
            
            readCountNew = readCountNew+1; % new contiguous portion means another count for newMatrix
            
            try
                % store observation in newMatrix and record start of new CpG observation
                newMatrix(readCountNew,observation) = dataMat(readNum,observation);
                CpGstart(readCountNew)              = observation; 
                
            catch % exceed newMatrix dimension, increase matrix size first and try again
                
                newMatrix = [newMatrix;-1*ones(2*origNumReads,size(dataMat,2))]; %#ok<AGROW>
                CpGstart  = [CpGstart;zeros(2*origNumReads,1)]; %#ok<AGROW>
                CpGend    = [CpGend;zeros(2*origNumReads,1)]; %#ok<AGROW>
                
                newMatrix(readCountNew,observation) = dataMat(readNum,observation);
                CpGstart(readCountNew)              = observation;
                
            end %end try catch
            
            if ~firstObservationOnRead
                CpGend(readCountNew-1) = prevObservation; % record the end of the previous read
            else
                firstObservationOnRead=0; % next iteration of loop will no longer be the first observation on the read
            end

            prevObservation = observation; % set up for next iteration of loop
            
        end %end test for contiguous portion
        
    end %end loop through observed CpGs in current read
    
    CpGend(readCountNew)=observation; % the last observation's index will be the end of the current read
        
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resize final output to only have relevant information (i.e., deal with
% the unused preallocated space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
newMatrix = newMatrix(1:readCountNew,:);
CpGstart  = CpGstart(1:readCountNew);
CpGend    = CpGend(1:readCountNew);

