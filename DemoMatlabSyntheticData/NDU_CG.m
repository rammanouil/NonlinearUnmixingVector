function [X_NDU, F_NDU, t_NDU, K_NDU] = NDU_CG(S,R,V,lambda,mu,par,str)


[L, N] = size(S);

% 1) sigma calculation
% ---------------------
if strcmp(str,'Separable') && par~=0 
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
if strcmp(str,'Separable') && par~=0
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
elseif strcmp(str,'Separable') && par==0
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

[X_NDU, F_NDU, t_NDU, K_NDU] = NDU_kernel_CG(S,R,lambda,mu,0.05,200,K,Eg);


