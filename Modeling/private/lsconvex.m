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
% lsconvex;
% check convexity
%

if nmin>1, 
  convex=0; 
else
  convex=1;   
  for i=2:s-1,
    f12=(flist(i)-flist(i-1))/(alist(i)-alist(i-1));
    f13=(flist(i)-flist(i+1))/(alist(i)-alist(i+1));
    f123=(f13-f12)/(alist(i+1)-alist(i-1));
    if f123<0, 
      if prt>1, disp('not convex'); end;
      convex=0;
      break; 
    end;
  end;
  if prt>1 & convex, disp('convex'); end;
end;
