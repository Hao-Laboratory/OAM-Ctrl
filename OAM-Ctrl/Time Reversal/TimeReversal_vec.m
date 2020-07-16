function [amp,phs,plr] = TimeReversal_vec(pupilDiaPixNum,wavelength,NA,n,Dipole)
%TIMEREVERSAL calculates the pupil function using time reversal.
%
% Ref. 
% "[1] Wan C, Yu Y, Zhan Q. Diffraction-limited near-spherical focal spot 
% with controllable arbitrary polarization using single objective lens[J]. 
% Opt Express, 2018, 26(21): 27109-27117."
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Jun.08, 2020

x = Dipole.Coord.X;  % coordinate
y = Dipole.Coord.Y;
z = Dipole.Coord.Z;
psiX = Dipole.Psi.X;
psiY = Dipole.Psi.Y;
psiZ = Dipole.Psi.Z;
aX = Dipole.A.X;
aY = Dipole.A.Y;
aZ = Dipole.A.Z;
dType = Dipole.Type;

psfNum = length(x);

[xx,yy] = meshgrid(linspace(-1,1,pupilDiaPixNum));
[phi,rho] = cart2pol(xx,yy);

theta = asin((NA/n).*rho);
k0 = 2*pi*n/wavelength;
kx = k0*sin(theta).*cos(phi);
ky = k0*sin(theta).*sin(phi);
kz = -k0*cos(theta);

M11 = 1+(cos(theta)-1).*cos(phi).^2;
M21 = (cos(theta)-1).*cos(phi).*sin(phi);
M31 = -sin(theta).*cos(phi);
M12 = (cos(theta)-1).*cos(phi).*sin(phi);
M22 = 1+(cos(theta)-1).*sin(phi).^2;
M32 = -sin(theta).*sin(phi);
M13 = sin(theta).*cos(phi);
M23 = sin(theta).*sin(phi);
M33 = cos(theta);

Dx2Theta = cos(theta).*cos(phi);
Dx2Phi = -sin(phi);
Dy2Theta = cos(theta).*sin(phi);
Dy2Phi = cos(phi);

Theta2Ex = cos(theta).*cos(phi);
Theta2Ey = cos(theta).*sin(phi);
Theta2Ez = sin(theta);
Phi2Ex = -sin(phi);
Phi2Ey = cos(phi);
Phi2Ez = 0;

Px = zeros(pupilDiaPixNum);
Py = Px;
Pz = Px;

for ii = 1:psfNum
    switch dType(ii)
        case 'e'
            Dz2Phi = 0;
            Dz2Theta = sin(theta);
        case 'm'
            Dz2Phi = sin(theta);
            Dz2Theta = 0;
    end
    Ex = (Dx2Theta.*Theta2Ex + Dx2Phi.*Phi2Ex).*aX(ii).*exp(-1i.*psiX{ii}) + ...
        (Dy2Theta.*Theta2Ex + Dy2Phi.*Phi2Ex).*aY(ii).*exp(-1i.*psiY{ii}) + ...
        (Dz2Theta.*Theta2Ex + Dz2Phi.*Phi2Ex).*aZ(ii).*exp(-1i.*psiZ{ii});
    
    Ey = (Dx2Theta.*Theta2Ey + Dx2Phi.*Phi2Ey).*aX(ii).*exp(-1i.*psiX{ii}) + ...
        (Dy2Theta.*Theta2Ey + Dy2Phi.*Phi2Ey).*aY(ii).*exp(-1i.*psiY{ii}) + ...
        (Dz2Theta.*Theta2Ey + Dz2Phi.*Phi2Ey).*aZ(ii).*exp(-1i.*psiZ{ii});
    
    Ez = (Dx2Theta.*Theta2Ez + Dx2Phi.*Phi2Ez).*aX(ii).*exp(-1i.*psiX{ii}) + ...
        (Dy2Theta.*Theta2Ez + Dy2Phi.*Phi2Ez).*aY(ii).*exp(-1i.*psiY{ii}) + ...
        (Dz2Theta.*Theta2Ez + Dz2Phi.*Phi2Ez).*aZ(ii).*exp(-1i.*psiZ{ii});
    
    
    Px = (M11.*Ex + M12.*Ey + M13.*Ez)./sqrt(cos(theta)).*...
        exp(-1i.*(kx.*x(ii) + ky.*y(ii) + kz.*z(ii))) + Px;
    
    Py = (M21.*Ex + M22.*Ey + M23.*Ez)./sqrt(cos(theta)).*...
        exp(-1i.*(kx.*x(ii) + ky.*y(ii) + kz.*z(ii))) + Py;
    
    % Pz == 0;
    Pz = (M31.*Ex + M32.*Ey + M33.*Ez)./sqrt(cos(theta)).*...
        exp(-1i.*(kx.*x(ii) + ky.*y(ii) + kz.*z(ii))) + Pz;
    
    textwaitbar(ii, psfNum, 'Time Reversaling');
end
amp = sqrt(abs(Px).^2 + abs(Py).^2);

phsX = angle(Px);
phsY = angle(Py);
delta_phi = phsY-phsX;

phs = mod(-(phsX),2*pi);

plrX = abs(Px)./amp;
plrY = abs(Py).*exp(-1i.*delta_phi)./amp;
plrX(amp==0) = 0;
plrY(amp==0) = 0;
plr(:,:,1) = plrX;
plr(:,:,2) = plrY;
plr(:,:,3) = zeros(pupilDiaPixNum);

% plr(:,:,1) = conj(Px);
% plr(:,:,2) = conj(Py);
% plr(:,:,3) = zeros(pupilDiaPixNum);

amp(rho>=1) = 0;
amp = amp./max(amp(:));
phs(rho>=1) = 0;
end