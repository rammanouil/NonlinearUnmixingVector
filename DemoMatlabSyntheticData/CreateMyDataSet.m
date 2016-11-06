function [L,X,S_lin,F0,F,S] = CreateMyDataSet(MM,R,N,P,coef,u,SNR)

[L,~]=size(R);

X = randg(coef,P, N);
X = X./(ones(P,1)*sum(X)); % abundances
S_lin = R*X; % Linear part in S

sigmaf = 1; 
widthf = round((6*sigmaf - 1)/2)-1;
h = fspecial('gaussian',[2*widthf+1 2*widthf+1],sigmaf);

sigmaf2 = 3.5; % 0.175*L; % 3.5; % 
h2 = fspecial('gaussian',[L 1],sigmaf2);
h2 = h2./max(h2(:));
h2 = repmat(h2,1,N);

if MM == 1 % LMM + bilinear contribution  
    F0 = u*S_lin.*S_lin; %
    F = F0; %
elseif MM == 2 % LMM + (bilinear + adjacency effect)  
    F0 = u*S_lin.*S_lin; % nonlinear contribution in a pixel
    F = conv2(padarray(F0,[widthf widthf],'replicate'), h,'valid');
elseif MM == 3 % LMM + (bilinear + adjacency effect + band selectivity) 
    F0 = u*S_lin.*S_lin; % nonlinear contribution in a pixel
    F1 = conv2(padarray(F0,[widthf widthf],'replicate'), h,'valid');
    F = F1.*h2;
else
    disp('error: model can be only assigned values 1 or 2');
    exit;
end

S = S_lin + F;
Pw_S = norm(S,'fro')^2/(N*L);
Pw_E = Pw_S/(10^(SNR/10));
sigma = sqrt(Pw_E);
S = S + sigma*randn(size(S));
