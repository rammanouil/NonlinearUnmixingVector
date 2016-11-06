
function [X, Z, F, t_ADMM, residues, K] = SuUNDU_Ki_kernel2_CG(S,R,lambda,mu,rho,nitermax,tol_pcg,nitermax_pcg,K,Eg,str2)

% 2) initialising
% ----------------
[L, N] = size(S);
P = size(R,2); 
A = [  speye(P); ones(1,P)];   % cst. var. % replace last row by zeros to remove sum-to-one
B = [ -speye(P); zeros(1,P)];  % cst. var. 
C = [0*speye(P,N); ones(1,N)]; % cst. var. % replace last row by zeros to remove sum-to-one 
invAtA = (-1/(P+1))*ones(P,P)+eye(P);
RinvAtARt = R*invAtA*R';

% if issparse(K)
%     kron_ISA = (1/rho)*kron(speye(N),R*invAtA*R'); % reserves less memory than full() version
%     Q = speye(L*N) + (1/lambda)*K + kron_ISA;
%     Q = (Q+Q')/2; 
% else
%     kron_ISA = (1/rho)*kron(eye(N),R*invAtA*R'); % reserves less memory than full() version
%     Q = eye(L*N) + (1/lambda)*K + kron_ISA;
%     Q = (Q+Q')/2;
% end
% clear kron_ISA; 
% X = zeros(P,N); % zeros(P,N); % For supervised & unsupervised case

[X, ~] = FCLS(S,R); 
Z = X;
Z_old = Z;
V = zeros(P+1,N);
niter = 1;
res = 1;
tol = sqrt(N)*1e-15;

residues = cell(1,2);
for i=1:2
    residues{1,i} = zeros(1,nitermax); % 1 primal residual, 2 dual
end

% ADMM algorithm 
% ---------------
t = clock;
w = zeros(L*N,1); 
if strcmp(str2,'UNDU')
    while (res > tol) && (niter <= nitermax)
        
        % {X,f} step
        % ----------
        p = vec(S + (1/rho)*R*invAtA*A'*(V+rho*B*Z-rho*C));
        % w = conjgrad(Q,p,w,1e-4,1000) ; % 5
        % [w, flag] = cgs(Q,p,1e-4,1000,[],[],w) ; % 5
        % [w, ~] = pcg(Q,p,tol_pcg,nitermax_pcg,[],[],w) ; % 5
        % w = Q\p;
        w = conjgrad_NDU(K,RinvAtARt,Eg,p,w,lambda,rho,tol_pcg,nitermax_pcg);
        W = reshape(w,L,N); % E = W ; E is equal to W
        X = (1/rho)*invAtA*(R'*W - A'*V - rho*A'*(B*Z-C));
        
        % Z step - In the unsupervised consists of the group lasso
        % --------------------------------------------------------
        for ii=1:P
            tmp = X(ii,:) + V(ii,:)/rho;
            tmp_plus = max(tmp,0);
            norm_tmp_plus = sqrt(sum(tmp_plus.^2));
            if norm_tmp_plus < (mu/rho)
                Z(ii,:) = zeros(size(tmp));
            else
                Z(ii,:) = (1 - (mu/rho)/norm_tmp_plus)*tmp_plus;
            end
        end
        
        % Update Lagrange multipliers
        % ----------------------------
        V = V + rho*(A*X + B*Z - C);
        
        % Stopping criteria
        % -----------------
        res_primal = norm(A*X+B*Z-C,'fro'); % rk
        res_dual = norm(A'*B*(Z_old - Z),'fro'); % sk
        res = res_primal+res_dual;
        Z_old = Z;
        
        % saving residuals
        % -----------------
        residues{1,1}(niter) = res_primal; % primal res. convergence
        residues{1,2}(niter) = res_dual;   % dual res. convergence
        
        % iter count
        % -----------
        niter = niter + 1;
        if mod(niter,20)
            disp(niter)
        end
    end
elseif strcmp(str2,'SuNDU') % ADMM 
    while (res > tol) && (niter <= nitermax)
        
        % {X,f} step
        % ----------
        p = vec(S + (1/rho)*R*invAtA*A'*(V+rho*B*Z-rho*C));
        % w = conjgrad(Q,p,w,1e-4,1000) ; % 5
        % [w, flag] = cgs(Q,p,1e-4,1000,[],[],w) ; % 5
        % [w, ~] = pcg(Q,p,tol_pcg,nitermax_pcg,[],[],w) ; % 5
        % w = Q\p;
        w = conjgrad_NDU(K,RinvAtARt,Eg,p,w,lambda,rho,tol_pcg,nitermax_pcg);
        W = reshape(w,L,N); % E = W ; E is equal to W
        X = (1/rho)*invAtA*(R'*W - A'*V - rho*A'*(B*Z-C));
        
        % Z step - In the supervised consists of projecting on R+
        % --------------------------------------------------------
        Z = max((rho/(rho+mu))*X+(1/(rho+mu))*V(1:end-1,:),0);
        
        % Update Lagrange multipliers
        % ----------------------------
        V = V + rho*(A*X + B*Z - C);
        
        % Stopping criteria
        % -----------------
        res_primal = norm(A*X+B*Z-C,'fro'); % rk
        res_dual = norm(A'*B*(Z_old - Z),'fro'); % sk
        res = res_primal+res_dual;
        Z_old = Z;
        
        % saving residuals
        % -----------------
        residues{1,1}(niter) = res_primal; % primal res. convergence
        residues{1,2}(niter) = res_dual;   % dual res. convergence
        
        % iter count
        % -----------
        niter = niter + 1;
        if mod(niter,100000000000)==0
            disp(niter)
        end
    end
elseif strcmp(str2,'GDNU') % QP 
    
    M = P ; 
    
    % Kv = eye(L*N)+ kron(K,Eg)/lambda + kron(R*R',eye(N))/mu;
    Kv = eye(L*N)+ K/lambda + kron(R*R',eye(N))/mu;
    A = kron([R sum(R,2)],speye(N))/mu;
    B = kron([speye(M) , ones(M,1); ones(1,M) , M]/mu,speye(N));
    
    Q = [Kv , A; A' , B];
    clear Kv A B;
    Q = (Q+Q')/2;
    % figure; spy(Q);
    
    p = [vec(S); zeros(M*N,1); ones(N,1)];
    
    tmp = (L*N+1):(L*N+M*N);
    Aineq = sparse(1:M*N,tmp,ones(1,M*N),M*N,N*(L+M+1),M*N);
    bineq = sparse(1:M*N,ones(M*N,1),zeros(M*N,1),M*N,1);
    
    % t=clock;
    % qpas: min 0.5 x'*Q*x + p'*x such that Ax < b
    [x] = qpas(full(Q),-p,full(-Aineq),full(bineq));
    % t = etime(clock,t);
    
    J = reshape(x(1:L*N),L,N);
    O = reshape(x(L*N+1:L*N+M*N),M,N);
    u = x(L*N+M*N+1:end);
    
    % f_est = kron(K,Eg)*vec(V')/lambda; % F_est = reshape(f_est,N,L)';
    % F = Eg*J*K/lambda; % K*V*Eg'/lambda;
    f = K*vec(J)/lambda; F = reshape(f,L,N);
    X = (R'*J+O+ones(M,1)*u')/mu;
end

t_ADMM = etime(clock,t);

if strcmp(str2,'GDNU')==0
    % F = K*w/lambda; 
    % F = reshape(F,L,N); % no need for this value at each iteration
    W = reshape(w,L,N); 
    F = (Eg*W*K)/lambda; 
end
