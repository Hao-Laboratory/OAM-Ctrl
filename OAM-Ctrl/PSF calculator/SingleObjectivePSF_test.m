% This is a test script for PSF calculation of single objective

clear; clc;
close all;

obj.NA = 1.4;  % numerical aperture (NA) of objective
obj.n = 1.518;  % refractive index of immersion medium

beam.wavelength = 1;  % wavelength of laser
beam.amp = ones(101);  % amplitude
beam.phs = vortexPhasePlate(101,2);  % phase
beam.plr = cylindVecPolarization(101,'l');  % polarization
beam.abr = zeros(101);  % aberration

% scope of calculation
scope.xs = linspace(-1,1,51);
scope.ys = linspace(-1,1,51);
scope.zs = 0;

% calculation
[Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
    101);

I = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;

figure;
pupilshowplr(beam.plr,15);
title('polarization');

figure;
subplot(2,2,1);
pupilshow(mod(beam.phs,2*pi));title('phase');
subplot(2,2,2);
imagesc(I);axis image xy off;colormap(gca,'hot');title('PSF');
subplot(2,2,3);
imagesc(abs(Ez).^2);axis image xy off;title('Iz');
subplot(2,2,4);
imagesc(mod(angle(Ez),2*pi));axis image xy off;title('phase');