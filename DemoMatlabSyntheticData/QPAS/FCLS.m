function [X_FCLS, t_FCLS] = FCLS(S,R)

H = R'*R; 
f = -R'*S;
P = size(R,2);
N = size(S,2);
l = zeros(P,1); 
A = ones(1,P);
b = 1;
X_FCLS = zeros(P,N);
t = clock;
for i=1:N
    % [x,err,lm] = qp(H,f,L,k,A,b,l,u,display);
    X_FCLS(:,i) = qpas(H,f(:,i),[],[],A,b,l);
end
t_FCLS = etime(clock,t);