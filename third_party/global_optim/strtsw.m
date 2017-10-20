%
% Permission was provided by the author, Arnold Neumaier, to 
% modify and distribute this MATLAB code with the informME 
% package. Contact the original author directly for use outside 
% this package. 
%
% Author's website:
% http://www.mat.univie.ac.at/~neum/
%
% Author's source websites:
% http://www.mat.univie.ac.at/~neum/software/mcs/
% http://www.mat.univie.ac.at/~neum/software/minq/
% http://www.mat.univie.ac.at/~neum/software/ls/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strtsw.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s = strtsw(smax,level,f)
% updates the record list for starting a new sweep and computes the
% lowest level containing non-split boxes

% Input:
% smax         depth of search
% level(1:nboxes)  levels of these boxes (level(j) = 0 if box j has 
%              already been split)
% f(1:nboxes)  their function values at the base vertices

% Output:
% s            lowest level with record(s) ~= 0

function s = strtsw(smax,level,f)
global nboxes record
% nboxes       number of boxes not in the `shopping basket'
% record(1:smax-1)  vector pointing to the best non-split box at each
%              level; record(s) = 0 if there is no non-split box at 
%              level s (record list)

% initialization
record = zeros(smax-1,1);
s = smax;
for j = 1:nboxes
  if level(j) > 0
    if level(j) < s
      s = level(j);
    end
    if ~record(level(j))  
      record(level(j)) = j;
    elseif f(j) < f(record(level(j)))
      record(level(j)) = j;
    end
  end
end
 
