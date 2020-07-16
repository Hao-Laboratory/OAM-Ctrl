function [amp,phs] = timeReversal(pupilDiaPixNum,wavelength,NA,n,coord,psi,A)
%TIMEREVERSAL calculates the pupil function using time reversal
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Apr.28, 2020

x = coord.X;
y = coord.Y;
z = coord.Z;
psfNum = length(x);

[xx,yy] = meshgrid(linspace(-1,1,pupilDiaPixNum));
[phi,rho] = cart2pol(xx,yy);

theta = asin((NA/n).*rho);
k0 = 2*pi*n/wavelength;
kx = k0*sin(theta).*cos(phi);
ky = k0*sin(theta).*sin(phi);
kz = -k0*cos(theta);

E = zeros(pupilDiaPixNum);
for ii = 1:psfNum
        E = A(ii)*exp(-1i*(kx.*(x(ii)) + ky.*(y(ii)) + kz.*(z(ii)) + psi{ii})).*sin(theta)./sqrt(cos(theta)) + E;
    textwaitbar(ii, psfNum, 'Time Reversaling');
end

P = conj(E);
amp = abs(P);
phs = mod(angle(P),2*pi);

amp(rho>1) = 0;
phs(rho>1) = 0;
end