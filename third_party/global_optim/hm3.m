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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hm3.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = hm3(x)
% Hartman3 function
function f = hm3(x)
if size(x,1) == 1
  x = x';
end
a = [3.0d0,  0.1d0,  3.0d0,  0.1d0;
     10.0d0, 10.0d0, 10.0d0, 10.0d0;
     30.0d0, 35.0d0, 30.0d0, 35.0d0];
p = [ 0.36890d0, 0.46990d0, 0.10910d0, 0.03815d0;
      0.11700d0, 0.43870d0, 0.87320d0, 0.57430d0;
      0.26730d0, 0.74700d0, 0.55470d0, 0.88280d0];
c = [1.0d0, 1.2d0, 3.0d0, 3.2d0];
for i=1:4
 d(i) = sum(a(:,i).*(x - p(:,i)).^2);
end
f = -sum(c.*exp(-d)); 
