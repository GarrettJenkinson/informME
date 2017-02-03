%
% Permission was provided by the author, Ali Mahammad-Djafari, 
% to modify and distribute this MATLAB code with the informME 
% package. Contact the original author directly for use outside 
% this package.
%
% Author's website:
% http://djafari.free.fr/index.htm
% 

function [lambda,p,entr] = maxent(mu,x,lambda0)

%
% ME_DENS2
% [LAMBDA,P,ENTR]=ME_DENS2(MU,X,LAMBDA0)
% This program calculates the Langrange Multipliers of the ME
% probability density functions p(x) from the knowledge of the
% N moment contstraints in the form:
% E{x^n}=mu(n) n=0:N with mu(0)=1.
%
% MU is a table containing the constraints MU(n),n=1:N.
% X is a table defining the range of the variation of x.
% LAMBDA0 is a table containing the first estimate of the LAMBDAs.
% (This argument is optional.)
% LAMBDA is a table containing the resulting Langrange parameters.
% P is a table containing the resulting pdf p(x).
% ENTR is a table containing the entropy values at each
% iteration.
%
% Author: A. Mohammad-Djafari
% Date : 10-01-1991
%

warning('off', 'MATLAB:illConditionedMatrix')
warning('off','MATLAB:nearlySingularMatrix')

maxIters = 5000;
entr = zeros(maxIters,1);

mu   = mu(:); 
mu   = [1;mu];
x    = x(:); 
lx   = length(x);
xmin = x(1); 
xmax = x(lx); 
dx   = x(2)-x(1);

if(nargin == 2)
    lambda = zeros(size(mu));
    lambda(1) = log(xmax-xmin);
else
    lambda = lambda0(:);
end
N = length(lambda);

M = 2*N-1;
fin = zeros(length(x),M);
fin(:,1) = ones(size(x));
for n = 2:M
    fin(:,n) = x.*fin(:,n-1);
end

iter = 0;
while 1
    iter = iter+1;
    p = exp(-(fin(:,1:N)*lambda));
    G = zeros(M,1);
    
    for n=1:M
        G(n) = dx*sum(fin(:,n).*p);
    end
    
    entr(iter) = lambda'*G(1:N);
    gnk = zeros(N,N);
    
    for i = 1:N
        gnk(:,i) = -G(i:N+i-1);
    end
    
    v      = mu-G(1:N);
    delta  = gnk\v; 
    lambda = lambda+delta;
    eps    = 1e-6;
    
    if(abs(delta./lambda)<eps), break, end
    
    if(iter>2)
        if(abs((entr(iter)-entr(iter-1))/entr(iter))<eps),break, end
    end
    
    if (iter>maxIters)
        disp('Warning: Maximum number of iterations reached'); break;
    end
    
end

p = exp(-(fin(:,1:N)*lambda));

entr = entr(1:iter);
warning('on', 'MATLAB:illConditionedMatrix')
warning('on','MATLAB:nearlySingularMatrix')

end
