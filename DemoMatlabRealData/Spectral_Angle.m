function A = Spectral_Angle(X1,X2)

% A : angles in degrees

corr = X1'*X2;
mod1 = sqrt(sum(X1.^2));
mod2 = sqrt(sum(X2.^2));

mod = mod1*mod2; 

A = acos(corr./(mod));
