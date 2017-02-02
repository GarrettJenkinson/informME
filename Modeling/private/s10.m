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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% s10.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s10(x)
% Shekel10 function
function f = s10(x)
a = [4, 1, 8, 6, 3, 2, 5, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 5, 1, 2, 3.6;
     4, 1, 8, 6, 3, 2, 3, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6];
c = [ 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5];
if size(x,1) == 1
 x = x';
end
for i=1:10
 b = (x - a(:,i)).^2;
 d(i) = sum(b);
end
f = -sum((c+d).^(-1));
