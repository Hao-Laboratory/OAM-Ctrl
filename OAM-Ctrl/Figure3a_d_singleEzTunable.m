% This script generates single longitudinal polarization optical
% vortex with tunable longitudinal strength in image space.
%
% The results are corresponding to Figure 3a-d in the paper.
%
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Dec.04, 2020

clear; clc;
close all;

%% variables
A = 0;  % ratio of radial polarization
l = 2;  % topological charges

%%
saveData = 0;

pupilDiaPixNum = 432;
pixelSize = 0.02;
quiverNumber = 20;

beam.wavelength = 1;
obj.NA = 0.95; obj.n = 1;
cmap_E = buildcmap('wbyr');
cmap_S = buildcmap('wkcr');
phi = vortexPhasePlate(pupilDiaPixNum,l);

Dipole.Coord.X = [0, 0];
Dipole.Coord.Y = [0, 0];
Dipole.Coord.Z = [0, 0];

Dipole.Type = ['e', 'm'];

Dipole.ThetaP = [0, 0];
Dipole.PhiP = [0, 0];

Dipole.Psi.X = {1*phi, 1*phi};
Dipole.Psi.Y = {1*phi, 1*phi};
Dipole.Psi.Z = {1*phi, 1*phi};

Dipole.A.X = [0, 0];
Dipole.A.Y = [0, 0];
Dipole.A.Z = [sqrt(A), sqrt(1-A)];

xHalfScope = 1.5;
beam.abr = zeros(pupilDiaPixNum);


%% time reversal
% [amp,phs,plr] = VectorialTimeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,Dipole);
[amp,phs,plr] = TimeReversal_vec(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,Dipole);

beam.amp = amp;
beam.phs = phs;
beam.plr = plr;

%% show pupil
figure;
subplot(131);
pupilshowplr(beam.plr,quiverNumber);
title('polarization');

subplot(132);
pupilshow(beam.amp);colormap(gca,'jet');
title('amplitude');

subplot(133);
pupilshow(beam.phs);
title('phase');

%% PSF calculation
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

% magnetic field strength and S
[Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum);
[Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
S = sqrt(Sx.^2+Sy.^2+Sz.^2);
Sxy_norm = sqrt(Sx.^2+Sy.^2)./max(S(:));

%% Plot
figure;
subplot(2,3,1);
imagesc(scope.xs,scope.xs,PSF_norm);
axis image xy off;colormap(gca,cmap_E);caxis([0,1]);title('{\itI_{xyz}} = 1');

subplot(2,3,2);
imagesc(scope.xs,scope.xs,Ir_norm);title(['{\itI_r} = ',num2str(round(max(Ir_norm(:)),2))]);
axis image xy off;colormap(gca,cmap_E);caxis([0,1]);

subplot(2,3,3);
imagesc(scope.xs,scope.xs,Iz_norm);title(['{\itI_z} = ',num2str(round(max(Iz_norm(:)),2))]);
axis image xy off;colormap(gca,cmap_E);caxis([0,1]);

subplot(2,3,4);
imagesc(mod(angle(Ez),2*pi));axis image xy off;title('\it\psi_z');

subplot(2,3,5);
imagesc(scope.xs,scope.ys,Sxy_norm);title('\itS_{xy}');
axis image xy off ;caxis([0 1]);
colormap(gca,cmap_S);
hold on;
myquiver2d(scope.xs,scope.ys,Sx,Sy,quiverNumber);

%% save data
if saveData == 1
    save('Data_singleEzTunalbe')
end