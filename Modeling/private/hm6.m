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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hm6.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = hm6(x)
% Hartman6 function

function f = hm6(x)
if size(x,1) == 1
 x = x';
end
a = [10.00,  0.05,  3.00, 17.00;
     3.00, 10.00,  3.50,  8.00;
     17.00, 17.00,  1.70,  0.05;
     3.50,  0.10, 10.00, 10.00;
     1.70,  8.00, 17.00,  0.10;
     8.00, 14.00,  8.00, 14.00];
p = [0.1312, 0.2329, 0.2348, 0.4047;
     0.1696, 0.4135, 0.1451, 0.8828;
     0.5569, 0.8307, 0.3522, 0.8732;
     0.0124, 0.3736, 0.2883, 0.5743;
     0.8283, 0.1004, 0.3047, 0.1091;
     0.5886, 0.9991, 0.6650, 0.0381];
c = [1.0, 1.2, 3.0, 3.2];
for i=1:4
 d(i) = sum(a(:,i).*(x - p(:,i)).^2);
end
f = -sum(c.*exp(-d)); 
