clear;
clc;
close all;

obj.NA = 1.4;
obj.n = 1.518;
beam.wavelength = 1;
beam.abr = zeros(101);
beam.amp = ones(101);

beam.phs = vortexPhasePlate(101,1);
beam.plr = CircularPolarization(101,'l');
scope.xs = linspace(-1,1,51);
scope.ys = linspace(-1,1,51);
scope.zs = 0;

[Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
    101,0,1,1);

I = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
figure;
subplot(1,3,1);
imagesc(I);axis image xy off;colormap(gca,'hot');
subplot(1,3,2);
imagesc(mod(beam.phs,2*pi));axis image xy off;
subplot(1,3,3);
imagesc(mod(angle(Ez),2*pi));axis image xy off;