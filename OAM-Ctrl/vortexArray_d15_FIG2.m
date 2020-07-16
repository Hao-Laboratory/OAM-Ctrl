% this script generates arbitrary longitudinal polarization optical
% vortices in random position with desired topological charges in image
% space by function modulation
%
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% May.05, 2020

clear;
close all;

%% common preset parameters
useIte = 0;
noIte = 0;

useAmp = 1;
saveData = 1;
ifSxy = 0;
xp = -7; yp = -13; zp = 0;

intensityError = 0.01;
stdThreshold = 0.0005;

memoryUnit = 500;
quiverScale = 10;
pupilDiaPixNum = 500;
pixelSize = 0.02;

if3D = 0;
showAllPsf = 1;

%% make folder
foldName = '12 points 3D anomaly intensity';

if saveData == 1
    if useIte == 1
        typeOption = 'ite';
    elseif noIte == 1
        typeOption = 'no ite';
    elseif noIte == 0
        typeOption = [];
    end
    
    switch useAmp
        case 1
            modOption = 'amp phs';
        case 0
            modOption = 'phs';
    end
    
    if if3D == 0
        folderName = ['.\Inverse Focusing\Time Reversal\demo\OAM ctrl\',...
            'pupil mod\data\arbitrary\',foldName,'\',modOption,'\',typeOption];
    elseif if3D == 1
        folderName = ['.\Inverse Focusing\Time Reversal\demo\OAM ctrl\',...
            'pupil mod\data\arbitrary\',foldName,'\',modOption,'\3D\',typeOption];
    end
    
    if showAllPsf == 1
        folderName = ['.\Inverse Focusing\Time Reversal\demo\OAM ctrl\',...
            'pupil mod\data\arbitrary\',foldName,'\',modOption,'\all Vortices'];
    end
    
    [status,msg] = mkdir(folderName);
    fileName = [folderName,'\Data.mat'];
end

%% the coordinates, amplitude, and phase of 13 points in two planes
xHalfScope = 15;
yHalfScope = 15;

% z = 0
coord1.X = [-7 11 0 -7 12];
coord1.Y = [11 11 0 -13 -13];
coord1.Z = zeros(1,5);

% z = 50
coord2.X = [-7 12 0 -12 11 4 -7];
coord2.Y = [12 11 1 -2 -4 -13 -11];
coord2.Z = 100*ones(1,7);

coord.X = [coord1.X, coord2.X];
coord.Y = [coord1.Y, coord2.Y];
coord.Z = [coord1.Z, coord2.Z];

l = 1;
phi = vortexPhasePlate(pupilDiaPixNum,l);

A0 = [1 0.2 0.4 1 1 1 1 0.8 0.5 0.7 0.8 1];
psi = {-4*phi, 3*phi, 5*phi, -2*phi, 2*phi, -5*phi, -4*phi, 4*phi, 6*phi,...
    -3*phi, 2*phi, 3*phi};

ringRadius = 2;

%% constant parameters
cmap_E = buildcmap('wbyr');
cmap_S = buildcmap('wkcr');

% objective lens
obj.NA = 0.95;
obj.n = 1;

% input beam
beam.wavelength = 1;
beam.abr = zeros(pupilDiaPixNum);
beam.plr = cylindVecPolarization(pupilDiaPixNum,0);

%% use iteration
if useIte == 1
    A_desired = A0;
    intensityThreshold = A_desired.*(1-intensityError);
    
    if saveData == 1
        gifName = [folderName,'\Animation.gif'];
    end
    
    %% PSF XY
    x = coord.X;
    y = coord.Y;
    z = coord.Z;
    
    psfPixNum = length(x);
    localMax = zeros(1,psfPixNum);
    STD = ones(1,1);
    f_iterationProcess = figure('outerposition',get(0,'screensize'));
    iteStep = 0;
    
    ns = length(-ringRadius:pixelSize:ringRadius);
    PSF = zeros(ns,ns,psfPixNum);
    psfMax = zeros(1,psfPixNum);
    xs = cell(1,ns); ys = cell(1,ns);
    while 1
        iteStep = iteStep +1;
        
        % time reversal
        [amp,phs] = timeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,...
            obj.n,coord,psi,A0);
        if useAmp == 1
            beam.amp = amp;
        elseif useAmp == 0
            beam.amp = ones(pupilDiaPixNum);
        end
        beam.phs = phs;
        
        % PSF calculate
        
        for ii = 1:psfPixNum
            scope.xs = x(ii)-ringRadius:pixelSize:x(ii)+ringRadius;
            scope.ys = y(ii)-ringRadius:pixelSize:y(ii)+ringRadius;
            scope.zs = z(ii);
            xs{ii} = scope.xs; ys{ii} = scope.ys;
            
            [xx,yy] = meshgrid(scope.xs,scope.ys);
            [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum,...
                1,memoryUnit,1);
            
            PSF(:,:,ii) = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
            psfMax(ii) = max(max(PSF(:,:,ii)));
        end
        PSF_norm = PSF./max(psfMax);
        
        for iii = 1:psfPixNum
            localMax(iii) = max(max(PSF_norm(:,:,iii)));
            
            if localMax(iii) < A_desired(iii)
                
                if std(localMax-A_desired) < 0.02
                    intensityFactor = 0.1;
                elseif std(localMax-A_desired) < 0.05
                    intensityFactor = 0.3;
                elseif std(localMax-A_desired) < 0.1
                    intensityFactor = 0.7;
                else
                    intensityFactor = 1;
                end
                
                A0(iii) =  A0(iii) + intensityFactor*(A_desired(iii)-...
                    localMax(iii));
            end
            
        end
        
        STD(iteStep) = std(localMax-A_desired);
        
        figure(f_iterationProcess);
        for jj = 1:psfPixNum
            subplot(4,4,jj);
            imagesc(xs{jj},ys{jj},PSF_norm(:,:,jj));
            title(['(',num2str(x(jj)),',',num2str(y(jj)),',',...
                num2str(z(jj)),')',' ',num2str(round(abs(A_desired(jj)-...
                localMax(jj)),2))]);
            axis image xy off;colormap(cmap_E);caxis([0,1]);
        end
        
        subplot(4,4,[13,14,15,16]);
        plot(iteStep,STD(iteStep),'*');
        hold on;grid on;
        xlabel('iterations');ylabel('STD');
        if saveData == 1
            drawnow;
            rgbImg = frame2im(getframe(f_iterationProcess));
            [indImg,colorMap] = rgb2ind(rgbImg,256);
            if iteStep == 1
                imwrite(indImg,colorMap,gifName,'gif', 'Loopcount',Inf,...
                    'DelayTime',0.8);
            else
                imwrite(indImg,colorMap,gifName,'gif','WriteMode',...
                    'append','DelayTime',0.8);
            end
        end
        
        % stop case
        if all(abs(1-localMax./A_desired) <= intensityError)
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
end

%% no iteration
if noIte == 1
    [amp,phs] = timeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,...
        coord,psi,A0);
    if useAmp == 1
        beam.amp = amp;
    elseif useAmp == 0
        beam.amp = ones(pupilDiaPixNum);
    end
    beam.phs = phs;
    
    load('vor_amp.mat');
    beam.amp = amp;
    load('vor_phs.mat');
    beam.phs = phs;
    
    scope.xs = -xHalfScope:pixelSize:xHalfScope;
    scope.ys = -yHalfScope:pixelSize:yHalfScope;
    
    for ii = 1:2
        if ii == 1
            if if3D == 0
                scope.zs = 0;
            elseif if3D == 1
                scope.zs = -3:pixelSize:3;
            end
            %% electric field
            [Ex_z0,Ey_z0,Ez_z0] = singleobjectivepsf(obj,beam,scope,...
                pupilDiaPixNum,1,memoryUnit,1);
            PSF_z0 = abs(Ex_z0).^2+abs(Ey_z0).^2+abs(Ez_z0).^2;
            PSF_z0_norm = PSF_z0./max(PSF_z0(:));
            
        elseif ii == 2
            if if3D == 0
                scope.zs = 100;
            elseif if3D == 1
                scope.zs = 100-3:pixelSize:100+3;
            end
            [Ex_z100,Ey_z100,Ez_z100] = singleobjectivepsf(obj,beam,scope,...
                pupilDiaPixNum,1,memoryUnit,1);
            PSF_z100 = abs(Ex_z100).^2+abs(Ey_z100).^2+abs(Ez_z100).^2;
            PSF_z100_norm = PSF_z100./max(PSF_z100(:));
        end
    end
end

%% show each PSF
if showAllPsf == 1
    load('vor_amp.mat');
    beam.amp = amp;
    load('vor_phs.mat');
    beam.phs = phs;
    
    N = length(coord.X);
    marks = {'A','B','C','D','E','F','G','H','I','J','K','L'};

    for ii = 1:N
        scope.xs = coord.X(ii)-2:pixelSize:coord.X(ii)+2;
        scope.ys = coord.Y(ii)-2:pixelSize:coord.Y(ii)+2;
        scope.zs = coord.Z(ii);

        %% electric field
        [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
            pupilDiaPixNum,1,memoryUnit,1);
        PSF(:,:,ii) = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
        Iz(:,:,ii) = abs(Ez).^2;
        
        %{
        %% magnetic field
        [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum,1,...
            memoryUnit,1);

        %% Poynting vector
        [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
        S = sqrt(Sx.^2+Sy.^2+Sz.^2);
        Sxy = sqrt(Sx.^2+Sy.^2);
        Sxy_norm = Sxy./max(S(:));
        Sxy_norm_max = max(Sxy_norm(:));
        Sxy_norm_min = min(Sxy_norm(:));
        
        f_phs_Ez = figure;
        imagesc(mod(angle(Ez),2*pi));axis image xy off;
        
        f_Sxy = figure;
        imagesc(scope.xs,scope.ys,Sxy_norm);axis image xy off;colormap(gca,cmap_S);
        caxis([0 1]);
%         colorbar('Ticks',[Sxy_norm_min,Sxy_norm_max],...
%             'TickLabels',{'0',num2str(round(Sxy_norm_max,2))});
        hold on;
        MyQuiver(scope.xs,scope.ys,Sx,Sy,quiverScale);
        
        mark = marks{ii};
        if saveData == 1
            figFormat = '-dpdf';
            print(f_phs_Ez,[folderName,'\phs Ez ',mark],figFormat,'-r300');
            print(f_Sxy,[folderName,'\Sxy ',mark],figFormat,'-r300');
        end
        %}
    end
    
    for jj = 1:N
        PSF_norm = PSF(:,:,jj)./max(PSF(:));
        Iz_norm = Iz(:,:,jj)./max(PSF(:));
        
        f_PSF = figure;
        imagesc(PSF_norm);axis image xy off;colormap(cmap_E);colorbar;caxis([0 1]);
        
        f_Ez = figure;
        imagesc(Iz_norm);axis image xy off;colormap(cmap_E);colorbar;caxis([0 1]);
        
        mark = marks{jj};
        if saveData == 1
            print(f_PSF,[folderName,'\PSF ',mark],figFormat,'-r300');
            print(f_Ez,[folderName,'\Ez ',mark],figFormat,'-r300');
        end
    end
    
end

%% Sxy
if ifSxy == 1
    scopeSize = 4;
    scope.xs = xp-scopeSize/2:pixelSize:xp+scopeSize/2;
    scope.ys = yp-scopeSize/2:pixelSize:yp+scopeSize/2;
    scope.zs = zp;
    
    [amp,phs] = timeReversal(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,...
        coord,psi,A0);
    if useAmp == 1
        beam.amp = amp;
    elseif useAmp == 0
        beam.amp = ones(pupilDiaPixNum);
    end
    beam.phs = phs;
    load('vor_amp.mat');
    beam.amp = amp;
    load('vor_phs.mat');
    beam.phs = phs;
    
    %% electric field
    [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum,1,...
        memoryUnit,1);
    
    %% magnetic field
    [Hx,Hy,Hz] = singleobjectivepsf_H(obj,beam,scope,pupilDiaPixNum,1,...
        memoryUnit,1);
    
    %% Poynting vector
    [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz);
    S = sqrt(Sx.^2+Sy.^2+Sz.^2);
    Sxy = sqrt(Sx.^2+Sy.^2);
    Sxy_norm = Sxy./max(S(:));
    Sxy_norm_max = max(Sxy_norm(:));
    Sxy_norm_min = min(Sxy_norm(:));
end

%% show amplitude and phase on pupil
f_pupil = figure;
subplot(1,2,1);
pupilshow(beam.amp);colormap(gca,'jet');
if useAmp == 1
    colorbar('Ticks',[min(beam.amp(:)),max(beam.amp(:))],...
        'TickLabels',{'0',num2str(round(max(beam.amp(:)),2))});
elseif useAmp == 0
    colorbar;
end

subplot(1,2,2);
pupilshow(phs);
colorbar('Ticks',[min(phs(:)),max(phs(:))],...
    'TickLabels',{'0',num2str(round(max(phs(:)),2))});

%% figure is no iteration
if noIte == 1
    if if3D == 0
        % z0
        % Ixyz
        Ixyz_z0 = PSF_z0_norm;
        
        f_Ixyz_z0 = figure;
        imagesc(scope.xs,scope.ys,Ixyz_z0);axis image xy;colormap(gca,cmap_E);
        colorbar('Ticks',[min(Ixyz_z0(:)),max(Ixyz_z0(:))],...
            'TickLabels',{'0',num2str(round(max(Ixyz_z0(:)),2))});
        set(gca,'XMinorTick','on','YMinorTick','on',...
            'XAxisLocation','origin','YAxisLocation','origin');
        
        % Iz
        Iz_z0 = abs(Ez_z0).^2./max(PSF_z0(:));
        
        f_Iz_z0 = figure;
        imagesc(scope.xs,scope.ys,Iz_z0);axis image xy off;colormap(gca,cmap_E);
        colorbar('Ticks',[min(Iz_z0(:)),max(Iz_z0(:))],...
            'TickLabels',{'0',num2str(round(max(Iz_z0(:)),2))});
        
        % phase of Ez
        phs_EZ_z0 = mod(angle(Ez_z0),2*pi);
        
        f_phsEz_z0 = figure;
        imagesc(scope.xs,scope.ys,phs_EZ_z0);axis image xy off;
        colorbar('Ticks',[min(phs_EZ_z0(:)),max(phs_EZ_z0(:))],...
            'TickLabels',{'0',num2str(round(max(phs_EZ_z0(:)),2))});
        
        % z100
        Ixyz_z100 = PSF_z100_norm;
        
        f_Ixyz_z100 = figure;
        imagesc(scope.xs,scope.ys,Ixyz_z100);axis image xy;colormap(gca,cmap_E);
        colorbar('Ticks',[min(Ixyz_z100(:)),max(Ixyz_z100(:))],...
            'TickLabels',{'0',num2str(round(max(Ixyz_z100(:)),2))});
        set(gca,'XMinorTick','on','YMinorTick','on',...
            'XAxisLocation','origin','YAxisLocation','origin');
        
        % Iz
        Iz_z100 = abs(Ez_z100).^2./max(PSF_z100(:));
        
        f_Iz_z100 = figure;
        imagesc(scope.xs,scope.ys,Iz_z100);axis image xy off;colormap(gca,cmap_E);
        colorbar('Ticks',[min(Iz_z100(:)),max(Iz_z100(:))],...
            'TickLabels',{'0',num2str(round(max(Iz_z100(:)),2))});
        
        % phase of Ez
        phs_EZ_z100 = mod(angle(Ez_z100),2*pi);
        
        f_phsEz_z100 = figure;
        imagesc(scope.xs,scope.ys,phs_EZ_z100);axis image xy off;
        colorbar('Ticks',[min(phs_EZ_z100(:)),max(phs_EZ_z100(:))],...
            'TickLabels',{'0',num2str(round(max(phs_EZ_z100(:)),2))});
    end
end

%% show Sxy
if ifSxy == 1
    % Sxy of focal field
    f_Sxy = figure('Color','w');
    imagesc(scope.xs,scope.ys,Sxy_norm);axis image xy off;colormap(gca,cmap_S);
    colorbar('Ticks',[Sxy_norm_min,Sxy_norm_max],...
        'TickLabels',{'0',num2str(round(Sxy_norm_max,2))});
    hold on;
    MyQuiver(scope.xs,scope.ys,Sx,Sy,quiverScale);
    set(gca,'BoxStyle','full');
end

%% save data
if saveData == 1
    save(fileName);
    figFormat = '-dpdf';
    
    print(f_pupil,[folderName,'\pupil'],figFormat,'-r300');
    
    if noIte == 1
        if if3D == 0
            print(f_Ixyz_z0,[folderName,'\Ixyz z0'],figFormat,'-r300');
            print(f_Iz_z0,[folderName,'\Iz z0'],figFormat,'-r300');
            print(f_phsEz_z0,[folderName,'\phs Ez z0'],figFormat,'-r300');
            
            print(f_Ixyz_z100,[folderName,'\Ixyz z100'],figFormat,'-r300');
            print(f_Iz_z100,[folderName,'\Iz z100'],figFormat,'-r300');
            print(f_phsEz_z100,[folderName,'\phs Ez z100'],figFormat,'-r300');
        end
    end
end