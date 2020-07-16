% this script generates bi-vortex structures in arbitrary position in image
% space controlled by pure phase modulation
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Apr.05, 2020

clear;
close all;

cmap_E = buildcmap('wbyr');
cmap_S = buildcmap('wkcr');
ifSaveData = 0;

ifXY = 1;
ifXZ = 0;

useAmp = 0; biAmp = 0;
if useAmp == 1
    modType = 'amp phs';
elseif useAmp == 0
    modType = 'phs';
end
memoryUnit = 500;
pupilDiaPixNum = 101;
pixelSize = 0.02;
quiverScale = 10;
xHalfScope = 2;
yHalfScope = 2;
l = 2; vortexDir = 'opposite';

x1 = 2; y1 = 0; z1 = -7;
x2 = -1; y2 = 0; z2 = 12;

beam.plr = cylindVecPolarization(pupilDiaPixNum,0);

%% optical parameters
obj.NA = 0.95;
obj.n = 1;

beam.wavelength = 1;
beam.amp = ones(pupilDiaPixNum);
beam.abr = zeros(pupilDiaPixNum);

%% phase generation
coord.X = [x1 x2];
coord.Y = [y1 y2];
coord.Z = [z1 z2];
phi = vortexPhasePlate(pupilDiaPixNum,l);

switch vortexDir
    case 'opposite'
        psi = {phi, -phi};
    case 'same'
        psi = {phi, phi};
end
A = [1 1];
[amp,phs] = timeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,coord,psi,A);

beam.phs = mod(phs,2*pi);

if useAmp == 1
    beam.amp = amp;
    if biAmp == 1
        beam.amp =  BinarizeAmplitude(amp);
    end
end

% show pupil
f_pupil = figure('Color','w');
subplot(1,2,1);
pupilshow(beam.amp);
colorbar;
colormap(gca,'jet');

subplot(1,2,2);
pupilshow(beam.phs);
colorbar('Ticks',[min(beam.phs(:)),max(beam.phs(:))],...
    'TickLabels',{'0','2\pi'});

%% PSF calculation
if ifXY == 1
    % z1
    scope.xs = x1-xHalfScope:pixelSize:x1+xHalfScope;
    scope.ys = y1-yHalfScope:pixelSize:y1+yHalfScope;
    scope.zs = z1;
    [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
        pupilDiaPixNum,1,memoryUnit,1);
    
    PSF = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2; PSF_max = max(PSF(:));
    PSF_norm = Normalization(PSF);
    
    Iz = abs(Ez).^2;
    Iz_norm = Iz./PSF_max; Iz_norm_max = max(Iz_norm(:)); Iz_norm_min = min(Iz_norm(:));
    phs_Ez = mod(angle(Ez),2*pi);
    
    [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum,1,memoryUnit,1);
    [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
    S = sqrt(Sx.^2+Sy.^2+Sz.^2);
    Sxy = sqrt(Sx.^2+Sy.^2);
    Sxy_norm = Sxy./max(S(:)); Sxy_norm_max = max(Sxy_norm(:)); Sxy_norm_min = min(Sxy_norm(:));
    
    f_PSF_z1_XY = figure;
    imagesc(scope.xs,scope.ys,PSF_norm);axis image xy off;colormap(gca,cmap_E);
    colorbar('Ticks',[0,1],...
        'TickLabels',{'0','1'});
    
    
    f_Ez_z1_XY = figure;
    imagesc(scope.xs,scope.ys,Iz_norm);axis image xy off;colormap(gca,cmap_E);
    caxis([0,1]);colorbar('Ticks',[0,1],'TickLabels',{'0','1'});
%     colorbar('Ticks',[Iz_norm_min,Iz_norm_max],...
%         'TickLabels',{'0',num2str(round(Iz_norm_max,2))});
    
    f_phs_z1_XY = figure;
    imagesc(scope.xs,scope.ys,phs_Ez);axis image xy off;
    
    f_Sxy_z1_XY = figure;
    imagesc(scope.xs,scope.ys,Sxy_norm);axis image xy off;colormap(gca,cmap_S);
    caxis([0,1]);colorbar('Ticks',[0,1],'TickLabels',{'0','1'});
%     colorbar('Ticks',[Sxy_norm_min,Sxy_norm_max],...
%         'TickLabels',{'0',num2str(round(Sxy_norm_max,2))});
    hold on;
    MyQuiver(scope.xs,scope.ys,Sx,Sy,quiverScale);
    
    % z2
    scope.xs = x2-xHalfScope:pixelSize:x2+xHalfScope;
    scope.ys = y2-yHalfScope:pixelSize:y2+yHalfScope;
    scope.zs = z2;
    [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
        pupilDiaPixNum,1,memoryUnit,1);
    
    PSF = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2; PSF_max = max(PSF(:));
    PSF_norm = Normalization(PSF);
    
    Iz = abs(Ez).^2;
    Iz_norm = Iz./PSF_max; Iz_norm_max = max(Iz_norm(:)); Iz_norm_min = min(Iz_norm(:));
    phs_Ez = mod(angle(Ez),2*pi);
    
    [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum,1,memoryUnit,1);
    [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
    S = sqrt(Sx.^2+Sy.^2+Sz.^2);
    Sxy = sqrt(Sx.^2+Sy.^2);
    Sxy_norm = Sxy./max(S(:)); Sxy_norm_max = max(Sxy_norm(:)); Sxy_norm_min = min(Sxy_norm(:));
    
    f_PSF_z2_XY = figure;
    imagesc(scope.xs,scope.ys,PSF_norm);axis image xy off;colormap(gca,cmap_E);
    
    f_Ez_z2_XY = figure;
    imagesc(scope.xs,scope.ys,Iz_norm);axis image xy off;colormap(gca,cmap_E);
    caxis([0,1]);
    
    f_phs_z2_XY = figure;
    imagesc(scope.xs,scope.ys,phs_Ez);axis image xy off;
    
    f_Sxy_z2_XY = figure;
    imagesc(scope.xs,scope.ys,Sxy_norm);axis image xy off;colormap(gca,cmap_S);
    caxis([0,1]);
    hold on;
    MyQuiver(scope.xs,scope.ys,Sx,Sy,quiverScale);
end

if ifXZ == 1
    % xz
    scope.xs = x2-xHalfScope-1:pixelSize:x1+xHalfScope+1;
    scope.ys = 0;
    scope.zs = z1-3:pixelSize:z2+3;
    
    [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
        pupilDiaPixNum,1,memoryUnit,1);
    
    PSF = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2; PSF_max = max(PSF(:));
    PSF_norm = Normalization(PSF);
    
    Ixy = abs(Ex).^2+abs(Ey).^2;
    Ixy_norm = Ixy./PSF_max; Ixy_norm_max = max(Ixy_norm(:)); Ixy_norm_min = min(Ixy_norm(:));
    Iz = abs(Ez).^2;
    Iz_norm = Iz./PSF_max; Iz_norm_max = max(Iz_norm(:)); Iz_norm_min = min(Iz_norm(:));
    
    f_PSF_XZ = figure;
    imagesc(scope.xs,scope.zs,PSF_norm);axis image xy off;colormap(gca,cmap_E);
    colorbar('Ticks',[0,1],'TickLabels',{'0','1'});

    f_Exy_XZ = figure;
    imagesc(scope.xs,scope.zs,Ixy_norm);axis image xy off;colormap(gca,cmap_E);
    colorbar('Ticks',[Ixy_norm_min,Ixy_norm_max],...
        'TickLabels',{'0',num2str(round(Ixy_norm_max,2))});
    
    f_Ez_XZ = figure;
    imagesc(scope.xs,scope.zs,Iz_norm);axis image xy;colormap(gca,cmap_E);
    caxis([0,1]);colorbar('Ticks',[0,1],'TickLabels',{'0','1'});
%     colorbar('Ticks',[Iz_norm_min,Iz_norm_max],...
%         'TickLabels',{'0',num2str(round(Iz_norm_max,2))});
end

%% save all data
if ifSaveData == 1
    folderName = ['.\Inverse Focusing\Time Reversal\demo\OAM ctrl\pupil mod\data\biVortex\',...
        num2str(vortexDir),'\',modType];
    [status,msg] = mkdir(folderName);
    
    
    figFormat = '-dpdf';
    
    print(f_pupil,[folderName,'\pupil'],figFormat,'-r300');
    
    print(f_PSF_z1_XY,[folderName,'\PSF_z1_XY'],figFormat,'-r300');
    print(f_Ez_z1_XY,[folderName,'\Ez_z1_XY'],figFormat,'-r300');
    print(f_phs_z1_XY,[folderName,'\phs_z1_XY'],figFormat,'-r300');
    print(f_Sxy_z1_XY,[folderName,'\Sxy_z1_XY'],figFormat,'-r300');
    
    print(f_PSF_z2_XY,[folderName,'\PSF_z2_XY'],figFormat,'-r300');
    print(f_Ez_z2_XY,[folderName,'\Ez_z2_XY'],figFormat,'-r300');
    print(f_phs_z2_XY,[folderName,'\phs_z2_XY'],figFormat,'-r300');
    print(f_Sxy_z2_XY,[folderName,'\Sxy_z2_XY'],figFormat,'-r300');
    
    print(f_PSF_XZ,[folderName,'\PSF_XZ'],figFormat,'-r300');
    print(f_Exy_XZ,[folderName,'\Exy_XZ'],figFormat,'-r300');
    print(f_Ez_XZ,[folderName,'\Ez_XZ'],figFormat,'-r300');
end