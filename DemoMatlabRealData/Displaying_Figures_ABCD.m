
clearvars ; close all; clc 

%% ============== Meris Image ======================= 
load('Golfe_Lion_ABCD2.mat')
M3d = ABCD_mer; clear ABCD_mer;

hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc((M3d(:,:,10))); colormap jet; colorbar 
AxisFont

%% ============== Class Image ======================= 
load('Golfe_Lion_ABCD_abd2');

class = unique(ABCD_class(:));
nb_class = length(class);

hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc((ABCD_class)); 
M = jet(nb_class); 
cbr = colormap(M); 
colormap(flipud(colormap));
colorbar('YTickLabel',{'121','142','222','243','313','323','332','512'}); AxisFont

%% ============== Endmembers ======================= 
clear all; 
load('Test_Final2.mat');

% E_det spectra
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
plot(R(:,[1 2 4 3]),'-*');
axis([1 L 0 0.4])
xlabel('Frequency band number');
ylabel('Reflectance');
legend('End. 1','End. 2','End. 3','End. 4');
AxisFont;

%% ============== Abd maps ======================= 
% FCLS 
for i=1:P
    tmp = reshape(X_FCLS(i,:),nl,nc);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 max(tmp(:))]); AxisFont
end

% ExtR
for i=1:P
    tmp = reshape(X_RExtR(i,:),nl,nc);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 max(tmp(:))]); AxisFont
end

% Khype (P)
for i=1:P
    tmp = reshape(X_Rkhype(i,:),nl,nc);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 1]); AxisFont
end

% Khype (G)
for i=1:P
    tmp = reshape(X_Rkhype(i,:),nl,nc);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 0.98],'YTickLabel',{'0','0.5','1'}); AxisFont
end

% GDNU  (P)
for i=1:P
    tmp = X3d2(:,:,i);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); %colormap jet; colorbar; AxisFont
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 0.98],'YTickLabel',{'0','0.5','1'}); AxisFont
end

% GDNU  (G)
for i=1:P
    tmp = X3d(:,:,i);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 0.98],'YTickLabel',{'0','0.5','1'}); AxisFont
end

%% ============== Nln maps ======================= 

% Khype (G) ********
for i=10
    tmp = reshape(F_Rkhype4(i,:),nl,nc);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); 
    M = jet; 
    cbr = colormap(M);  % colorbar  
    colorbar('YTick',[0 0.05 0.1]); AxisFont
end

% GDNU  (P) ********
for i=10
    tmp = F3d2(:,:,i);
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(tmp); %colormap jet; colorbar; AxisFont
    M = jet; 
    cbr = colormap(M); 
    %colorbar 
    colorbar('YTick',[0 0.05 0.1]); AxisFont
end

%% ============== Nln maps ======================= 
tmp = reshape(F_SuNDU2',nl,nc,L);
tmp_sub = tmp(171:180,231:240,:);
tmp_sub_ = reshape(tmp_sub,length(tmp_sub(:))/L,L)';
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(tmp_sub_); 
colorbar; AxisFont 
caxis manual ; caxis([0 0.12]); colormap jet; colorbar; AxisFont
colorbar('YTick',[0.01 0.06 0.1]); AxisFont

tmp = reshape(F_Rkhype4',nl,nc,L);
tmp_sub = tmp(171:180,231:240,:);
tmp_sub_ = reshape(tmp_sub,length(tmp_sub(:))/L,L)';
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(tmp_sub_); 
colorbar; AxisFont 
caxis manual ; caxis([0 0.12]); colormap jet; colorbar; AxisFont
colorbar('YTick',[0.01 0.06 0.1]); AxisFont

%% ============== Cstm Abd maps ======================= 
% Abundance maps deduced from the CLC classification maps 
clear all; 
load('Custumized_abd_ABCD'); 
for i=1:3
    hFig = figure;
    set(hFig, 'Position', [1 1000 300 200])
    imagesc(abd_map3(:,:,i)); 
    M = gray; 
    cbr = colormap(M); 
    colorbar('YTick',[0 0.5 0.98],'YTickLabel',{'0','0.5','1'}); AxisFont
end
