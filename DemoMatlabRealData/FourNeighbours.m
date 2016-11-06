function [S1, S2, S3, S4] = FourNeighbours(M3d,h,w,i,j)

[nl,nc,L] = size(M3d);

S1 = zeros(L,h*w); 
S2 = zeros(L,h*w); 
S3 = zeros(L,h*w); 
S4 = zeros(L,h*w); 

xi = (i-1)*h+1:i*h; 
yj = (j-1)*w+1:j*w;

cnt = 0;

for yp = yj
    for xp = xi
        
        cnt = cnt + 1;
        
        % left neighbour
        xl = xp;
        yl = max(yp-1,1);
        % right neighbour
        xr = xp;
        yr = min(yp+1,nc);
        % upper
        xu = max(xp-1,1);
        yu = yp;
        % lower
        xlo = min(xp+1,nl);
        ylo = yp;
        
        S1(:,cnt) = squeeze(M3d(xl,yl,:));
        S2(:,cnt) = squeeze(M3d(xr,yr,:));
        S3(:,cnt) = squeeze(M3d(xu,yu,:));
        S4(:,cnt) = squeeze(M3d(xlo,ylo,:));
        
    end
end








