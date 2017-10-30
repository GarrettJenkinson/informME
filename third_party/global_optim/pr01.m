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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pr01.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function summe=pr01(name,x)
% prints a (0,1) profile of x and returns the number of nonzeros
%
function summe=pr01(name,x)

% check for row or column vector
[m,n]=size(x);row=(m==1);col=(n==1);n=max(m,n); 
if ~(row|col), error('x must be a vector'); end;

text=[name,': '];summe=0;
for k=1:n, 
  if x(k), 
    text=[text,'1'];
    summe=summe+1;
  else 
    text=[text,'0']; 
  end;
end;
disp([text,'   ',num2str(summe),' nonzeros']);
