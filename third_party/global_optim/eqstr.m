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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eqstr.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function equal=eqstr(str1,str2)
% checks equality of two strings
%
function equal=eqstr(str1,str2)

equal=0;
if size(str1,2)~=size(str2,2), return; end;
equal=min(str1==str2);
