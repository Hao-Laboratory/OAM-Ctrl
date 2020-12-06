% This script generates arbitrary longitudinal polarization optical
% vortices in image space by pupil function modulation.
%
% The results are corresponding to Figure 2 and Figure S2 in the paper.
%
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Dec.04, 2020

clear; clc;
close all;

%% common preset parameters
if3D = 0;  % show 3D intensity distribution (Figure 2c)
saveData = 0;  % save data

intensityError = 0.01;  % intensity deviation between result and expectation 
stdThreshold = 0.0005;  % break when the standard deviation is stable

quiverNumber = 20;  % the number of arrows for energy flow illustration
pupilDiaPixNum = 432;  % resolution of the pupil function
pixelSize = 0.02;  % pixel size of image space

%% coordinates, amplitude, and phase of 12 vortices in two planes
xHalfScope = 15;  % range of image space
yHalfScope = 15;

marks = {'A','B','C','D','E','F','G','H','I','J','K','L'};  % marks of optical vortices

% desired positions
coord1.X = [-7 11 0 -7 12];
coord1.Y = [11 11 0 -13 -13];
coord1.Z = zeros(1,5);  % z = 0

coord2.X = [-7 12 0 -12 11 4 -7];
coord2.Y = [12 11 1 -2 -4 -13 -11];
coord2.Z = 100*ones(1,7);  % z = 100

coord.X = [coord1.X, coord2.X];
coord.Y = [coord1.Y, coord2.Y];
coord.Z = [coord1.Z, coord2.Z];

% topological charges
l = 1;
phi = vortexPhasePlate(pupilDiaPixNum,l);

I_desired = [1 0.2 0.4 1 1 1 1 0.8 0.5 0.7 0.8 1];
psi = {-4*phi, 3*phi, 5*phi, -2*phi, 2*phi, -5*phi, -4*phi, 4*phi, 6*phi,...
    -3*phi, 2*phi, 3*phi};
I0 = I_desired;

%% constant parameters
cmap_E = buildcmap('wbyr');
cmap_S = buildcmap('wkcr');

% objective lens
obj.NA = 1.4;
obj.n = 1.518;

% input beam
beam.wavelength = 1;  % wavelength
beam.abr = zeros(pupilDiaPixNum);  % aberration
beam.plr = cylindVecPolarization(pupilDiaPixNum,0);  % polarization

%% optimization process
x = coord.X;
y = coord.Y;
z = coord.Z;

psfPixNum = length(x);
localMax = zeros(1,psfPixNum);
STD = ones(1,1);
f_iterationProcess = figure('WindowState','maximized');
iteStep = 0;

% ROI of each vortex calculation
roi = 2;

% memory preallocation
ns = length(-roi:pixelSize:roi);
PSF = zeros(ns,ns,psfPixNum);
psfMax = zeros(1,psfPixNum);
Iz = PSF;
err = psfMax;
xs = cell(1,12); ys = xs;
Exx = cell(1,12); Eyy = Exx; Ezz = Exx;
while 1
    iteStep = iteStep + 1;
    
    % time reversal
    [amp,phs] = timeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,...
        obj.n,coord,psi,sqrt(I0));
    beam.amp = amp;
    beam.phs = phs;
    
    % PSF calculation
    for ii = 1:psfPixNum
        scopeROI.xs = x(ii)-roi:pixelSize:x(ii)+roi;
        scopeROI.ys = y(ii)-roi:pixelSize:y(ii)+roi;
        scopeROI.zs = z(ii);
        xs{ii} = scopeROI.xs; ys{ii} = scopeROI.ys;
        
        [xx,yy] = meshgrid(scopeROI.xs,scopeROI.ys);
        [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scopeROI,pupilDiaPixNum);
        Exx{ii} = Ex;
        Eyy{ii} = Ey;
        Ezz{ii} = Ez;
        PSF(:,:,ii) = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
        Iz(:,:,ii) = abs(Ez).^2;
        psfMax(ii) = max(max(PSF(:,:,ii)));
    end
    
    PSF_norm = PSF./max(psfMax);
    Iz_norm = Iz./max(psfMax);
    
    % optimize the amplitude of dipoles
    for iii = 1:psfPixNum
        localMax(iii) = max(max(PSF_norm(:,:,iii)));
        if localMax(iii) < I_desired(iii)
            err(iii) = I_desired(iii)-localMax(iii);
            intensityFactor = 1;
            I0(iii) =  I0(iii) + intensityFactor*(err(iii));
        end
    end
    
    STD(iteStep) = std(err);
    
    figure(f_iterationProcess);
    for ii = 1:psfPixNum
        subplot(4,4,ii);
        imagesc(xs{ii},ys{ii},PSF_norm(:,:,ii));
        title([marks{ii},'(',num2str(x(ii)),',',num2str(y(ii)),',',...
            num2str(z(ii)),')','  error: ',num2str(round(err(ii),3))]);
        axis image xy off;colormap(cmap_E);caxis([0,1]);
    end
    
    subplot(4,4,[13,14,15,16]);
    plot(1:iteStep,STD(1:iteStep));
    axis tight;
    grid on;
    xlabel('iterations');ylabel(['STD: ',num2str(round(STD(iteStep),3))]);
    drawnow;
    
    % stop case
    if all(abs(1-localMax./I_desired) <= intensityError)
        disp('Intensity is OK! Time Reversal Finished!');
        break;
    end
    
    if iteStep > 10
        if (std(STD(end-10:end)) < stdThreshold)
            disp('Standard Deviation is OK! Time Reversal Finished!');
            break;
        end
    end
end

%% show pupil function
f_pupil = figure;
subplot(1,2,1);
pupilshow(beam.amp);
colormap(gca,'jet');
colorbar;
title('amplitude');

subplot(1,2,2);
pupilshow(beam.phs);
colorbar;
title('phase');

%% overall PSF calculation
scope.xs = -xHalfScope:pixelSize:xHalfScope;
scope.ys = -yHalfScope:pixelSize:yHalfScope;

for ii = 1:2
    if ii == 1
        if if3D == 0
            scope.zs = 0;
        elseif if3D == 1
            scope.zs = -3:pixelSize:3;
        end
        [Ex_z0,Ey_z0,Ez_z0] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum);
        PSF_z0 = abs(Ex_z0).^2+abs(Ey_z0).^2+abs(Ez_z0).^2;
        Ixyz_z0 = PSF_z0./max(PSF_z0(:));
        Iz_z0 = abs(Ez_z0).^2./max(PSF_z0(:));
        phs_EZ_z0 = mod(angle(Ez_z0),2*pi);
        
    elseif ii == 2
        if if3D == 0
            scope.zs = 100;
        elseif if3D == 1
            scope.zs = 100-3:pixelSize:100+3;
        end
        [Ex_z100,Ey_z100,Ez_z100] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum);
        PSF_z100 = abs(Ex_z100).^2+abs(Ey_z100).^2+abs(Ez_z100).^2;
        Ixyz_z100 = PSF_z100./max(PSF_z100(:));
        Iz_z100 = abs(Ez_z100).^2./max(PSF_z100(:));
        phs_EZ_z100 = mod(angle(Ez_z100),2*pi);
    end
end

% plot
if if3D == 0
    fig_overallIntensity = figure('WindowState','maximized');
    subplot(2,3,1);
    imagesc(scope.xs,scope.ys,Ixyz_z0);axis image xy;colormap(gca,cmap_E);
    caxis([0 1]);colorbar;ylabel('{\itz} = 0');title('\itI');
    
    subplot(2,3,2);
    imagesc(scope.xs,scope.ys,Iz_z0);axis image xy off;colormap(gca,cmap_E);
    caxis([0 1]);colorbar;title('\itI_z');
    
    subplot(2,3,3);
    imagesc(scope.xs,scope.ys,phs_EZ_z0);axis image xy off;
    colorbar;title('\it\psi_z');
    
    subplot(2,3,4);
    imagesc(scope.xs,scope.ys,Ixyz_z100);axis image xy;colormap(gca,cmap_E);
    caxis([0 1]);colorbar;ylabel('{\itz} = 100');
    
    subplot(2,3,5);
    imagesc(scope.xs,scope.ys,Iz_z100);axis image xy off;colormap(gca,cmap_E);
    caxis([0 1]);colorbar;
    
    subplot(2,3,6);
    imagesc(scope.xs,scope.ys,phs_EZ_z100);axis image xy off;
    colorbar;
    
elseif if3D == 1
    fig_overallIntensity = figure;
    h = volrender(scope.xs,scope.ys,80-3:pixelSize:80+3,Iz_z0,cmap_E,linspace(0,0.1,101));
    hold on;
    volrender(scope.xs,scope.ys,scope.zs,Iz_z100,cmap_E,linspace(0,0.1,101));
    Ax = h.parent;
    Ax.XLabel.String = '';
    Ax.YLabel.String = '';
    Ax.ZLabel.String = '';
    Ax.ZTick = [80 100];
    Ax.XTickLabel = {'','',''};
    Ax.YTickLabel = {'','',''};
    Ax.ZTickLabel = {'','',''};
    Ax.LineWidth = 3;
end

%% single vortex calculation
N = length(coord.X);

% memory preallocation
ls = length(-roi:pixelSize:roi);
phsEz = zeros(ls,ls,N);
fig_singleVortex = figure('WindowState','maximized');
for ii = 1:N
    Ex = Exx{ii};
    Ey = Eyy{ii};
    Ez = Ezz{ii};
    phs_Ez = mod(angle(Ez),2*pi);
    
    I = PSF_norm(:,:,ii);
    Iz = Iz_norm(:,:,ii);
    
    scope.xs = coord.X(ii)-roi:pixelSize:coord.X(ii)+roi;
    scope.ys = coord.Y(ii)-roi:pixelSize:coord.Y(ii)+roi;
    scope.zs = coord.Z(ii);
    
    % magnetic field
    [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum);
    
    % Poynting vector
    [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
    S = sqrt(Sx.^2+Sy.^2+Sz.^2);
    Sxy = sqrt(Sx.^2+Sy.^2);
    Sxy_norm = Sxy./max(S(:));
    Sxy_norm_max = max(Sxy_norm(:));
    Sxy_norm_min = min(Sxy_norm(:));
    
    mark = marks{ii};
    
    subplot(4,12,ii);
    imagesc(I);axis image xy off;colormap(gca,cmap_E);caxis([0 1]);
    title(mark);
    
    subplot(4,12,12+ii);
    imagesc(Iz);axis image xy off;colormap(gca,cmap_E);caxis([0 1]);
    
    subplot(4,12,24+ii);
    imagesc(phs_Ez);axis image xy off;
    
    subplot(4,12,36+ii);
    imagesc(scope.xs,scope.ys,Sxy_norm);axis image xy off;colormap(gca,cmap_S);
    caxis([0 1]);
    
    hold on;
    myquiver2d(scope.xs,scope.ys,Sx,Sy,quiverNumber);
end

%% save data
if saveData == 1
    save('Data_ArbitraryVortices.mat');
end