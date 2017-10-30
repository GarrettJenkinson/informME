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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quadmin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = quadmin(a,b,d,x0)
% computes the minimum x of the quadratic polynomial
% p(x) = d(1) + d(2)(x - x0(1)) + d(3)(x - x0(1))(x - x0(2))
% in the interval [a,b]
% Uses the following function/m-file:
% quadpol.m

function x = quadmin(a,b,d,x0)
if d(3) == 0
 if d(2) > 0
  x = a;
 else
  x = b;
 end
elseif d(3) > 0
 x1 = 0.5*(x0(1) + x0(2)) - 0.5*d(2)/d(3);
 if a <= x1 & x1 <= b
  x = x1;
 elseif quadpol(a,d,x0) < quadpol(b,d,x0)
  x = a;
 else
  x = b;
 end
else
 if quadpol(a,d,x0) < quadpol(b,d,x0)
  x = a;
 else
  x = b;
 end
end
