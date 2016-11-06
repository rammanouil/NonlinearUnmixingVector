function [avg_angle, cnt] = Compute_avg_angle(A,B)

N = size(A,2); 
cnt = 0; 
tmp = 0; 

for i=1:N
    x = Spectral_Angle(A(:,i),B(:,i)); % angle [0 pi]  
    if ~isnan(x)
        tmp = tmp + (x);
        cnt = cnt+1; 
    end
end

avg_angle = tmp/cnt; % [0 pi]
