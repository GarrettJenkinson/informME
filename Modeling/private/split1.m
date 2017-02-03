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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function z = split1(x1,x2,f1,f2)
% Input: points x1 and x2, x1 < x2, and corresponding function values f1
%        and f2
% splits the interval [x1,x2] according to the golden section rule
% the part containing the better point gets the larger fraction of the 
% interval

function z = split1(x1,x2,f1,f2)
if f1 <= f2
  z = x1 + 0.5*(-1 + sqrt(5))*(x2 - x1);
else
  z = x1 + 0.5*(3 - sqrt(5))*(x2 - x1);
end
