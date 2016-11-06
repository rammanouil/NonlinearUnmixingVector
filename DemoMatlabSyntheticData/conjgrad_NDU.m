
function [x] = conjgrad_NDU(K,D,E,b,x,lambda,rho,tol,nitermax)

% K,D,lambda,rho,p,zeros(L*N,1),1e-6,100
% here A = I + (1/lambda)KoE + (1/rho)I_NoD
% A*x has three terms : x + vec + vec 

L = size(D,1); 
N = size(K,1); 
% r = b - A*x; % residu  
X = reshape(x,L,N); 
r = b - ( x + (1/lambda)*vec(E*X*K) + (1/rho)*vec(D*X) ); 
% gradient descent direction = residu
p = r; 
rsold = r'*r;

for i=1:nitermax
    
    % Ap = A*p;
    P = reshape(p,L,N); 
    Ap = p + (1/lambda)*vec(E*P*K) + (1/rho)*vec(D*P); 
    
    alpha=rsold/(p'*Ap); % optimal descent step 

    x=x+alpha*p; % descent 
    
    r=r-alpha*Ap; % update residu 
    
    if sqrt(rsold)<tol
        break;
    end
    
    rsnew = (r'*r);
    
    beta = rsnew/rsold; 
    
    p=r+beta*p; % update descent direction 
    
    rsold=rsnew;
    
end 

end

% r=b-A*x; % residu  
% p=r; % gradient descent direction = residu
% rsold=r'*r;
%  
% for i=1:nitermax
%     Ap = A*p; 
%     alpha=rsold/(p'*Ap); % optimal descent step 
%     x=x+alpha*p; % descent 
%     r=r-alpha*Ap; % update residu 
%     if sqrt(rsold)<tol
%         break;
%     end
%     beta = (r'*r)/rsold; 
%     p=r+beta*p; % update descent direction 
%     rsold=rsnew;
% end
