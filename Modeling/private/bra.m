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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bra.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = bra(x,y)
% Branin's function

function f = bra(x,y)
a=1;
b=5.1/(4*pi*pi);
c=5/pi;
d=6;
h=10;
ff=1/(8*pi);
if nargin == 1
  x1 = x(1);
  x2 = x(2);
else
  x1 = x;
  x2 = y;
end
f=a.*(x2-b.*x1.^2+c.*x1-d).^2+h.*(1-ff).*cos(x1)+h;
