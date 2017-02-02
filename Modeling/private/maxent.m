% License permission provided by the author to Garrett Jenkinson 
% to modify and distribute with informME package. Contact
% the original author directly for use outside of informME.
%
% Author's website:
% http://djafari.free.fr/index.htm
% 
function [lambda,p,entr]=maxent(mu,x,lambda0)
%
%ME_DENS2
% [LAMBDA,P,ENTR]=ME_DENS2(MU,X,LAMBDA0)
% This program calculates the Lagrange Multipliers of the ME
% probability density functions p(x) from the knowledge of the
% N moment contstraints in the form:
% E{x^n}=mu(n) n=0:N with mu(0)=1.
%
% MU is a table containing the constraints MU(n),n=1:N.
% X is a table defining the range of the variation of x.
% LAMBDA0 is a table containing the first estimate of the LAMBDAs.
% (This argument is optional.)
% LAMBDA is a table containing the resulting Lagrange parameters.
% P is a table containing the resulting pdf p(x).
% ENTR is a table containing the entropy values at each
% iteration.
%
% Author: A. Mohammad-Djafari
% Date : 10-01-1991
%
% Modified by Garrett Jenkinson on 02-03-1015 with modifications commented: 
% '% ADDED BY GARRETT JENKINSON' or '% REMOVED BY GARRETT JENKINSON'
% To turn off warnings, uncomment the next lines % ADDED BY GARRETT JENKINSON
warning('off', 'MATLAB:illConditionedMatrix')  % ADDED BY GARRETT JENKINSON 
warning('off','MATLAB:nearlySingularMatrix')   % ADDED BY GARRETT JENKINSON 

maxIters = 5000; % ADDED BY GARRETT
entr = zeros(maxIters,1); % ADDED BY GARRETT JENKINSON

mu=mu(:); mu=[1;mu]; % add mu(0)=1
x=x(:); lx=length(x); % x axis
xmin=x(1); xmax=x(lx); dx=x(2)-x(1);
%
if(nargin == 2) % initialize LAMBDA
    lambda=zeros(size(mu)); % This produces a uniform
    lambda(1)=log(xmax-xmin); % distribution.
else
    lambda=lambda0(:);
end
N=length(lambda);
%
M=2*N-1; % Calcul de fin(x)=x.^n
fin=zeros(length(x),M); %
fin(:,1)=ones(size(x)); % fi0(x)=1
for n=2:M
    fin(:,n)=x.*fin(:,n-1);
end
%
iter=0;
while 1 % start iterations
    iter=iter+1;
    %
    p=exp(-(fin(:,1:N)*lambda)); % Calculate p(x)
    %plot(x,p); % plot it
    %
    G=zeros(M,1); % Calculate Gn
    for n=1:M
        G(n)=dx*sum(fin(:,n).*p);
    end
    %
    entr(iter)=lambda'*G(1:N); % Calculate the entropy value
    %
    gnk=zeros(N,N); % Calculate gnk
    for i=1:N % Matrix G is a Hankel matrix
        gnk(:,i)=-G(i:N+i-1);
    end
    %
    v=mu-G(1:N); % Calculate v
    delta=gnk\v; % Calculate delta
    lambda=lambda+delta; % Calculate lambda
    eps=1e-6; % Stopping rules
    if(abs(delta./lambda)<eps), break, end
    if(iter>2)
        if(abs((entr(iter)-entr(iter-1))/entr(iter))<eps),break, end
    end
    
    
    if (iter>maxIters) % ADDED BY GARRETT JENKINSON 
        disp('Warning: Maximum number of iterations reached'); break;% ADDED BY GARRETT JENKINSON
    end                % ADDED BY GARRETT JENKINSON
    
end
%
p=exp(-(fin(:,1:N)*lambda)); % Calculate the final p(x)

%entr=entr(:); % REMOVED BY GARRETT JENKINSON
entr=entr(1:iter); % ADDED BY GARRETT JENKINSON
warning('on', 'MATLAB:illConditionedMatrix')  % ADDED BY GARRETT JENKINSON
warning('on','MATLAB:nearlySingularMatrix')   % ADDED BY GARRETT JENKINSON  
end
