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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% neighbor.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2] = neighbor(x,delta,u,v)
% computes 'neighbors' x1 and x2 of x needed for making triple search 
% and building a local quadratic model such that x(i), x1(i), x2(i) are
% pairwise distinct for i = 1,...,n
% Input:
% x  	points for which 'neighbors' should be found
% delta	distance of the neighbors
% [u,v]	original box for global optimization

function [x1,x2] = neighbor(x,delta,u,v)
i1 = find(x == u);
i2 = find(x == v);	
x1 = max(u,x-delta);
x2 = min(x+delta,v);
x1(i1) = x(i1) + 2*delta(i1);
x2(i2) = x(i2) - 2*delta(i2);

