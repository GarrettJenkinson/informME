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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vert1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,x1,x2,f1,f2] = vert1(j,z,f,x1,x2,f1,f2)
% subprogram used by vertex.m

function [x,x1,x2,f1,f2] = vert1(j,z,f,x1,x2,f1,f2)
if j == 1
  j1 = 2;
else
  j1 = 1;
end
x = z(j1);
if x1 == Inf
  x1 = z(j);
  f1 = f1 + f(j);
elseif x2 == Inf & x1 ~= z(j)
  x2 = z(j);
  f2 = f2 + f(j);
end
