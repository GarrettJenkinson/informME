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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% range.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function range(u,v,x,p)
% computes the range [amin,amax] for making a line search along x + ap
% when x + ap is restricted to [u,v]
%
function [amin,amax] = range(u,v,x,p)
amin = -Inf; amax = Inf;
for i=1:length(x)
  if p(i)>0, 
    amin=max(amin,(u(i)-x(i))/p(i));
    amax=min(amax,(v(i)-x(i))/p(i));
  elseif p(i)<0, 
    amin=max(amin,(v(i)-x(i))/p(i));
    amax=min(amax,(u(i)-x(i))/p(i));
  end
end  
