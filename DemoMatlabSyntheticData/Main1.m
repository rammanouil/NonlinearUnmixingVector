
clearvars ; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Create Data set               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MM = 3; % mixing model 1 2 or 3
SNR = 30; % signal to noise ratio 
N = 100; % number of pixels 
stp = 20; % downsampling number of spectral bands 
coef = 1; % coeficient of beta distribution
u = 0.2; % attenuation coeficient of nonlinear term 
load('endmembers.mat'); % Dictionary of available endmembers
R = M(1:stp:end,[2 6 7]); 
% R = M(1:stp:end,[4 6 7 8]);
% R = M(1:stp:end,[4 6 7 8 9]);
P = size(R,2); % number of endmembers 
[L,X,S_lin,F0,F,S] = CreateMyDataSet(MM,R,N,P,coef,u,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Ext end. method                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RExt = zeros(L,P+(P^2-P)/2);
RExt(:,1:P) = R;
cnt = P+1;
for i= 1:P
    for j=i+1:P
        RExt(:,cnt) = R(:,i).*R(:,j);
        cnt = cnt+1;
    end
end
[X_Ext, t_Ext] = FCLS_no1(S,RExt);
F_Ext = RExt(:,P+1:end)*X_Ext(P+1:end,:);
[RMSE_X_Ext, std_X_Ext] = ErrComput(X_Ext(1:P,:),X);
[RMSE_F_Ext, std_F_Ext] = ErrComput(F_Ext,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    khype (P)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0.01;
lambda = 1;
tmp = pdist2(R,R);
par = 0; % 0: Polynomial kernel & 1: Gaussian kernel  
sigma = par*max(tmp(:));
t = clock; 
[X_khype, F_khype, K_Khype] = khype3(sigma,lambda,mu,N,S,R);
t_khype = etime(clock,t);
[RMSE_X_khype, std_X_khype] = ErrComput(X_khype,X);
[RMSE_F_khype, std_F_khype] = ErrComput(F_khype,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    NDU (Tr+P)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 10;
mu = 0.1;
str = 'Transformable'; % Separable or Transformable
par = 0; % 0: Polynomial kernel & 1: Gaussian kernel  
V = cell(1,N); % Neighbours 
S1 = circshift(S,[0 1]);
S2 = circshift(S,[0 -1]);
for i=1:N
    V{1,i} = [S(:,i) S1(:,i) S2(:,i)];
end
[X_NDU, F_NDU, t_NDU, K_NDU] = NDU(S,R,V,lambda,mu,par,str);
[RMSE_X_NDU, std_X_NDU] = ErrComput(X_NDU,X);
[RMSE_F_NDU, std_F_NDU] = ErrComput(F_NDU,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    NDU (Sp+P)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.01;
mu = 0.0001;
str = 'Separable'; % Separable or Transformable
par = 0; % 0: Polynomial kernel & 1: Gaussian kernel  
V = cell(1,N); % Neighbours 
S1 = circshift(S,[0 1]);
S2 = circshift(S,[0 -1]);
for i=1:N
    V{1,i} = [S(:,i) S1(:,i) S2(:,i)];
end
[X_NDU2, F_NDU2, t_NDU2, K_NDU2] = NDU(S,R,V,lambda,mu,par,str);
[RMSE_X_NDU2, std_X_NDU2] = ErrComput(X_NDU2,X);
[RMSE_F_NDU2, std_F_NDU2] = ErrComput(F_NDU2,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Figures                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(strcat('Ext: ',num2str(RMSE_X_Ext),' ; ',num2str(RMSE_F_Ext))); 
disp(strcat('khype: ',num2str(RMSE_X_khype),' ; ',num2str(RMSE_F_khype))); 
disp(strcat('NDU: ',num2str(RMSE_X_NDU),' ; ',num2str(RMSE_F_NDU))); 
disp(strcat('NDU2: ',num2str(RMSE_X_NDU2),' ; ',num2str(RMSE_F_NDU2))); 

F1 = F;
F2 = F_Ext;
F3 = F_khype;
F4 = F_NDU;
F5 = F_NDU2;

tst = [F1 F2 F3 F4];
minn = min(tst(:)); 
maxx = max(tst(:)); 
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F1); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F2); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F3); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F4); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F5); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont


