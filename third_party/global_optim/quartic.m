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
% function f=quartic(a,x);
% evaluates f(x)=a(1)x^4+a(2)x^3+a(3)x^2+a(4)x+a(5)
% simultaneously at a vector x
%
function f=quartic(a,x);

f=(((a(1)*x+a(2)).*x+a(3)).*x+a(4)).*x+a(5); 

