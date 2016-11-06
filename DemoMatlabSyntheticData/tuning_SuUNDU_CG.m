
function [res, res_bst, par_bst, X_RSuNDU, Z_RSuNDU, F_RSuNDU, t_RSuNDU, residues, K ] = tuning_SuUNDU_CG(lambda_,mu_,par, ...
    S,R,V,str,X,F,nitermax,str2)

n1 = length(lambda_);
n2 = length(mu_);

res = cell(1,3);
res{1,1} = zeros(n1,n2); % RMSE_X_Rkhype
res{1,2} = zeros(n1,n2); % RMSE_F_Rkhype
res{1,3} = zeros(n1,n2); % SA_Rkhype
[L, N] = size(S);

% 1) sigma calculation
% ---------------------
if strcmp(str,'GraphReg') && par~=0 
    sigma = zeros(N,N); 
    for i=1:N
        for j=1:N
            v_i = V{1,i}; % setting data at band i
            v_j = V{1,j}; % setting data at band i
            K = pdist2(v_i(:)',v_j(:)');
            sigma(i,j) = max(K(:));
        end
    end
    par = par*max(sigma(:)); % par becomes the value of sigma to use next 
end

% 2) Kernel calculation 
% ---------------------
if strcmp(str,'GraphReg') && par~=0
    % E
    Wg = 1*circshift(eye(L),[0 1]); Wg(end,1)=0; Wg = (Wg+Wg');
    Dg = diag(sum(Wg,2)+1);
    tmp = Dg-Wg; 
    Eg = pinv(tmp);
    % K
    K = zeros(N,N); % K = zeros(L*N,L*N); tic % 54 sec and 584820000  double Vs 0.414943 seconds and 6908408 double
    for i=1:N
        for j=1:N
            v_i = V{1,i}; % setting data at band i
            v_j = V{1,j}; % setting data at band i
            Kij = generate_kernel_2(v_i(:),v_j(:),par);
            %Kij = Kij;
            K(i,j) = Kij;
        end
    end
elseif strcmp(str,'GraphReg') && par==0
    % E
    Wg = 1*circshift(eye(L),[0 1]); Wg(end,1)=0; Wg = (Wg+Wg');
    Dg = diag(sum(Wg,2)+1);
    tmp = Dg-Wg; 
    Eg = pinv(tmp);
    % K
    K = zeros(N,N); % K = zeros(L*N,L*N); tic % 54 sec and 584820000  double Vs 0.414943 seconds and 6908408 double
    for i=1:N
        for j=1:N
            v_i = V{1,i}; % setting data at band i
            v_j = V{1,j}; % setting data at band i
            Kij = (v_i(:)'*v_j(:)).^2; % generate_kernel_2(v_i(:),v_j(:),par);
            %Kij = Kij;
            K(i,j) = Kij;
        end
    end
end

K = K./max(K(:));

for i_n1=1:n1
    for j_n2=1:n2
        tic 
        lambda = lambda_(i_n1); % non linear contribution
        mu = mu_(j_n2);  % linear contribution
        %[X_RSuNDU, Z_RSuNDU, F_RSuNDU, t_RSuNDU, residues, K] = SuUNDU_Ki_kernel2(S,R,lambda,mu,0.05,nitermax,...
        %   1e-4,1000,K,str2);
        [X_RSuNDU, Z_RSuNDU, F_RSuNDU, t_RSuNDU, residues, K] = SuUNDU_Ki_kernel2_CG(S,R,lambda,mu,0.05,nitermax,...
            1e-4,1000,K,Eg,str2);
        toc
        % Il n'y a pas de mauvaise graines ni de mauvais cultivateurs mais
        res{1,1}(i_n1,j_n2) = Compute_RMSE(X,X_RSuNDU);
        res{1,2}(i_n1,j_n2) = Compute_RMSE(F,F_RSuNDU);
        res{1,3}(i_n1,j_n2) = Compute_avg_angle(F,F_RSuNDU);
        
    end
end


[min1, idx1] = min(res{1,1}); 
[~, idx2] = min(min1); 
idx1 = idx1(idx2); 

if (n1~=1) && (n2~=1)
    lambda = lambda_(idx1);
    mu = mu_(idx2);
    [X_RSuNDU, Z_RSuNDU, F_RSuNDU, t_RSuNDU, residues, K] = SuUNDU_Ki_kernel2_CG(S,R,lambda,mu,0.05,nitermax,...
        1e-4,1000,K,Eg,str2);
end

res_bst = zeros(1,3);         
res_bst(1) = Compute_RMSE(X_RSuNDU,X);
res_bst(2) = Compute_RMSE(F_RSuNDU,F);
res_bst(3) = Compute_avg_angle(F_RSuNDU,F);
if strcmp(str2,'UNDU')
    res_bst(3) = sum(mean(Z_RSuNDU,2)>0.01);
end

par_bst = [lambda, mu]; 


