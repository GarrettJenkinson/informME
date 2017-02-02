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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shu.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = shu(x,y)
% Shubert's function
function f = shu(x,y)
sum1 = 0; sum2 = 0;
if nargin == 1
 x1 = x(1);
 x2 = x(2);
else
 x1 = x;
 x2 = y;
end
for i = 1:5
 sum1 = sum1 + i.*cos((i+1).*x1+i);
 sum2 = sum2 + i.*cos((i+1).*x2+i);
end
f = sum1.*sum2;
