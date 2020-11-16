function [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,...
    pupilDiaPixNum)
%SINGLEOBJECTIVEOSF_H calculates the magnetic field strength 
% 'H' in the vicinity of the focal spot
% 
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Apr.23, 2020

beamH = beam;
beamH.plr(:,:,1) = -beam.plr(:,:,2);  % Hx = -Ey;
beamH.plr(:,:,2) = beam.plr(:,:,1);  % Hy = Ex;
beamH.plr(:,:,3) = beam.plr(:,:,3);  % Hz = Ez;

[Hx,Hy,Hz] = singleobjectivepsf(obj,beamH,scope,...
    pupilDiaPixNum);

Epsilon = 8.854e-12;  % permittivity in vacuum
Mu = 4*pi*1e-7;  % permeability in vacuum
prefix = sqrt(Epsilon/Mu);

Hx = prefix*Hx;
Hy = prefix*Hy;
Hz = prefix*Hz;
end