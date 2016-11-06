

clear all ; 
load('modelP3Smoothx0Smoothf2SNR40SuNDUfullS.mat');
F00 = F;
F_RSuNDU_full00 = F_GDNU2;
tst = [F00 F_RSuNDU_full00];
minn = min(tst(:)); 
maxx = max(tst(:)); 
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F00); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont
hFig = figure;
set(hFig, 'Position', [1 1000 300 200])
imagesc(F_RSuNDU_full00); colorbar; AxisFont % caxis manual ; caxis([minn-eps maxx+eps]); colorbar; AxisFont




