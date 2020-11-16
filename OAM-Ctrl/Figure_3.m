% This script is used to optimize the strength of the longitudinal
% component of single optical vortex
%
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Nov.16, 2020

clear; clc;
close all;

pupilDiaPixNum = 101;
pixelSize = 0.1;

Iz_desired = 0.8;
l = -1;

beam.wavelength = 1;
obj.NA = 0.95; obj.n = 1;
cmap_E = buildcmap('wbyr');
cmap_S = buildcmap('wkcr');
phi = vortexPhasePlate(pupilDiaPixNum,l);

Dipole.Coord.X = [0, 0];
Dipole.Coord.Y = [0, 0];
Dipole.Coord.Z = [0, 0];

Dipole.Type = ['e', 'm'];

Dipole.Psi.X = {0, 0};
Dipole.Psi.Y = {0, 0};
Dipole.Psi.Z = {1*phi, 1*phi};

Dipole.A.X = [0 0];
Dipole.A.Y = [0 0];
Dipole.A.Z = [sqrt(Iz_desired), sqrt(1-Iz_desired)];

xHalfScope = 1.5;
beam.abr = zeros(pupilDiaPixNum);

%% calculation
[amp,phs,plr] = TimeReversal_vec(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,Dipole);

beam.amp = amp;
beam.phs = phs;
beam.plr = plr;

scope.xs = -xHalfScope:pixelSize:xHalfScope;
scope.ys = -xHalfScope:pixelSize:xHalfScope;
scope.zs = 0;

[Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum);

PSF = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
PSF_norm = PSF./max(PSF(:));

Ir = abs(Ex).^2+abs(Ey).^2;
Ir_norm = Ir./max(PSF(:));

Iz = abs(Ez).^2;
Iz_norm = Iz./max(PSF(:));

f_PSF_Ez = figure;
subplot(1,3,1);
imagesc(scope.xs,scope.xs,PSF_norm);
axis image xy off;colormap(cmap_E);caxis([0,1]);

subplot(1,3,2);
imagesc(scope.xs,scope.xs,Ir_norm);title(num2str(round(max(Ir_norm(:)),2)));
axis image xy off;colormap(cmap_E);caxis([0,1]);

subplot(1,3,3);
imagesc(scope.xs,scope.xs,Iz_norm);title(num2str(round(max(Iz_norm(:)),2)));
axis image xy off;colormap(cmap_E);caxis([0,1]);

%% show pupil
f_pupil = figure;
subplot(131);
pupilshowplr(beam.plr,50);
subplot(132);
pupilshow(beam.amp);colorbar;colormap(gca,'jet');
subplot(133);
pupilshow(beam.phs);colorbar;

%% magnetic field strength and S
[Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum);
[Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
S = sqrt(Sx.^2+Sy.^2+Sz.^2);
Sxy_norm = sqrt(Sx.^2+Sy.^2)./max(S(:));


%% plot
f_phs_Ez = figure;
imagesc(mod(angle(Ez),2*pi));axis image xy off;

f_Sxy = figure;
imagesc(scope.xs,scope.ys,Sxy_norm);
axis image xy off ;caxis([0 1]);
colormap(gca,cmap_S);
hold on;
MyQuiver(scope.xs,scope.ys,Sx,Sy,10);