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
% lssplit;		
% find split of interval [i,i+1] towards better point
% no lssort at the end!
%

if flist(i)<flist(i+1), fac=short;
elseif flist(i)>flist(i+1), fac=1-short;
else fac=0.5;
end;
alp=alist(i)+fac*(alist(i+1)-alist(i));
if prt>2, disp(['split at ',num2str(alp)]); end;


