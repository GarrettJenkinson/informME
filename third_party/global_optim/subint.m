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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2] = subint(x,y)
% computes for real x and real or infinite y two points x1 and x2 in
% [min(x,y),max(x,y)] that are neither too close nor too far away from x

function [x1,x2] = subint(x,y)
x2 = y;
f = 1000;
if f*abs(x) < 1
 if abs(y) > f, x2 = sign(y); end
else
 if abs(y) > f*abs(x), x2 = 10*sign(y)*abs(x); end
end
x1 = x + (x2 - x)/10;
