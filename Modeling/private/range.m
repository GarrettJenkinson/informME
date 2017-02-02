% License permission provided by the author to Garrett Jenkinson 
% to modify and distribute with informME package. Contact
% the original author directly for use outside of informME.
%
% Author website:
% http://www.mat.univie.ac.at/~neum/
% 
% Author source websites:
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
