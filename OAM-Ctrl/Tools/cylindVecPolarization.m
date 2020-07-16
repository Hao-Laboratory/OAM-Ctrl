function plr = cylindVecPolarization(pupilDiaPixNum,theta)
%VECPOLARIZATION generates cylinderical vector beams, theta is the angle 
% between the polarization vector and the radial vector.
%
% Especially, 
% theta = 0 or pi: radial polarization
% theta = pi/2 or -pi/2: azimuthal polarization
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Jun.09, 2020

x = linspace(-1,1,pupilDiaPixNum);
[xx, yy] = meshgrid(x);
[phi, ~] = cart2pol(xx,yy);

plr(:,:,1) = cos(phi+theta);
plr(:,:,2) = sin(phi+theta);
plr(:,:,3) = zeros(pupilDiaPixNum);
end