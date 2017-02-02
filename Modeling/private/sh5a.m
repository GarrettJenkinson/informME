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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sh5.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = sh5(x)
% Shekel5 function
function f = sh5(x)
a = [4.0d0, 1.0d0, 8.0d0, 6.0d0, 3.0d0;
     4.0d0, 1.0d0, 8.0d0, 6.0d0, 7.0d0;
     4.0d0, 1.0d0, 8.0d0, 6.0d0, 3.0d0;
     4.0d0, 1.0d0, 8.0d0, 6.0d0, 7.0d0];
c = [0.1d0, 0.2d0, 0.2d0, 0.4d0, 0.4d0];
for i=1:5
b = (x - a(:,i)).^2;
 d(i) = sum(b);
end
f = -sum((c+d).^(-1));
