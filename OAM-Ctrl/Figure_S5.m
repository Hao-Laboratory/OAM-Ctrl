% this script is used to optimize the longitudinal component of every
% optical vortex in the focused field
%
% Author: Xin Liu
% Email: liuxin2018@zju.edu.cn
% Nov.16, 2020

clear; clc
close all;
pupilDiaPixNum = 101;
pixelSize = 0.1;

stdThreshold = 0.0002;
intensityError = 0.01;

beam.wavelength = 1;
obj.NA = 0.95; obj.n = 1;
cmap_E = buildcmap('wbyr');
phi = vortexPhasePlate(pupilDiaPixNum,1);

d = 30;
coord.X = [-d 0 d -d 0 d -d  0  d];
coord.Y = [ d d d  0 0 0 -d -d -d];
N = length(coord.X);
coord.Z = zeros(1,N);

Dipole.Coord.X = [coord.X, coord.X];
Dipole.Coord.Y = [coord.Y, coord.Y];
Dipole.Coord.Z = [coord.Z, coord.Z];

Dipole.Type = [repmat('e',1,N),repmat('m',1,N),];

Dipole.Psi.X = mat2cell(zeros(1,2*N),1,ones(1,2*N));
Dipole.Psi.Y = mat2cell(zeros(1,2*N),1,ones(1,2*N));
Dipole.Psi.Z = {1*phi, -1*phi, 2*phi, 1*phi, -1*phi, 2*phi, 2*phi, 2*phi, 2*phi,...
    1*phi, -1*phi, 2*phi, 1*phi, -1*phi, 2*phi, 2*phi, 2*phi, 2*phi};

Dipole.A.X = zeros(1,2*N);
Dipole.A.Y = zeros(1,2*N);

A1 = sqrt([1, 1, 1, 0, 0, 0, 0.5, 0.6, 0.7]);

Dipole.A.Z = [A1, sqrt(1-A1.^2)];
ringRadius = 2;
beam.abr = zeros(pupilDiaPixNum);

%% use iteration
Ez_desired = (Dipole.A.Z(1:N)).^2;
A_desired = ones(1,N);
eff_ind = find( Dipole.A.Z(1:N)~=1 & Dipole.A.Z(1:N)~=0 );

%% PSF XY
x = coord.X;
y = coord.Y;
z = coord.Z;

psfPixNum = length(x);
localMax_I = zeros(1,psfPixNum);
localMax_Ez = zeros(1,psfPixNum);

STD_I = ones(1,1);
STD_Ez = ones(1,1);

f_iterationProcess = figure;
iteStep = 0;

ns = length(-ringRadius:pixelSize:ringRadius);
PSF = zeros(ns,ns,psfPixNum);
Ez_stack = zeros(ns,ns,psfPixNum);

psfMax = zeros(1,psfPixNum);
EzMax = zeros(1,psfPixNum);
Ez_norm = zeros(ns,ns,psfPixNum);
xs = cell(1,ns); ys = cell(1,ns);

[xxi,yyi] = meshgrid(-ringRadius:pixelSize:ringRadius);
[~,rho] = cart2pol(xxi,yyi);
while 1
    iteStep = iteStep +1;
    
    % time reversal
    [amp,phs,plr] = TimeReversal_vec(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,Dipole);
    
    beam.amp = amp;
    beam.phs = phs;
    beam.plr = plr;
    
    % PSF calculate
    for ii = 1:psfPixNum
        scope.xs = x(ii)-ringRadius:pixelSize:x(ii)+ringRadius;
        scope.ys = y(ii)-ringRadius:pixelSize:y(ii)+ringRadius;
        scope.zs = z(ii);
        xs{ii} = scope.xs; ys{ii} = scope.ys;
        
        [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum);
        
        PSF(:,:,ii) = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
        Ez_stack(:,:,ii) = abs(Ez).^2;
        psfMax(ii) = max(max(PSF(:,:,ii)));
        
        Ez_norm(:,:,ii) = Ez_stack(:,:,ii)./max(psfMax(ii));
        EzMax(ii) = max(max(Ez_norm(:,:,ii)));
    end
    
    % intensity optimization
    PSF_norm = PSF./max(psfMax);
    
    for iii = 1:psfPixNum
        localMax_I(iii) = max(max(PSF_norm(:,:,iii)));
        
        devStd = std(localMax_I-A_desired);
        devI = A_desired(iii)-localMax_I(iii);

            if devStd < 0.01
                intensityFactor_I = 0.1;
            elseif devStd < 0.05
                intensityFactor_I = 0.5;
            else
                intensityFactor_I = 1;
            end
            
            if Ez_desired(iii) == 1
                Dipole.A.Z(iii) =  Dipole.A.Z(iii) + ...
                    intensityFactor_I*sign(devI)*sqrt(abs(devI));
            elseif Ez_desired(iii) == 0
                Dipole.A.Z(iii+N) =  Dipole.A.Z(iii+N) + ...
                    intensityFactor_I*sign(devI)*sqrt(abs(devI));
            else
                Dipole.A.Z(iii) =  Dipole.A.Z(iii) + ...
                    0.5*intensityFactor_I*sign(devI)*sqrt(abs(devI));
                Dipole.A.Z(iii+N) =  Dipole.A.Z(iii+N) + ...
                    0.5*intensityFactor_I*sign(devI)*sqrt(abs(devI));
            end
    end
    
    % Ez optimization
    for iii = 1:length(eff_ind)
        localMax_Ez(eff_ind(iii)) = EzMax(eff_ind(iii));
        
        devStd = std(localMax_Ez(eff_ind)-Ez_desired(eff_ind));
        devE = Ez_desired(eff_ind(iii)) - localMax_Ez(eff_ind(iii));
            if devStd < 0.01
                intensityFactor_Ez = 0.1;
            elseif devStd < 0.05
                intensityFactor_Ez = 0.5;
            else
                intensityFactor_Ez = 1;
            end
            Dipole.A.Z(eff_ind(iii)) =  Dipole.A.Z(eff_ind(iii)) + ...
                intensityFactor_Ez*sign(devE)*sqrt(abs(devE));
            Dipole.A.Z(eff_ind(iii)+N) =  Dipole.A.Z(eff_ind(iii)+N) - ...
                intensityFactor_Ez*sign(devE)*sqrt(abs(devE));
    end
    
    STD_I(iteStep) = std(localMax_I-A_desired);
    STD_Ez(iteStep) = std(localMax_Ez(eff_ind)-Ez_desired(eff_ind));
    
    figure(f_iterationProcess);
    for jj = 1:psfPixNum
        subplot(4,N,jj);
        imagesc(xs{jj},ys{jj},PSF_norm(:,:,jj));
        title(['(',num2str(x(jj)),',',num2str(y(jj)),',',...
            num2str(z(jj)),')',' ',num2str(round(abs(A_desired(jj)-...
            localMax_I(jj)),3))]);
        axis image xy off;colormap(cmap_E);caxis([0,1]);
        
        subplot(4,N,N+jj);
        imagesc(xs{jj},ys{jj},Ez_norm(:,:,jj));
        title(['(',num2str(x(jj)),',',num2str(y(jj)),',',...
            num2str(z(jj)),')',' ',num2str(round(abs(Ez_desired(jj)-...
            localMax_Ez(jj)),3))]);
        axis image xy off;colormap(cmap_E);caxis([0,1]);
    end
    
    subplot(4,N,3*N:1:4*N);
    if iteStep > 1
        plot(iteStep-1:iteStep,[STD_I(iteStep-1),STD_I(iteStep)],'r*');
        hold on;grid on;
        plot(iteStep-1:iteStep,[STD_Ez(iteStep-1),STD_Ez(iteStep)],'k+');
        xlim([1 iteStep]);
        xlabel('iterations');ylabel('STD');
        set(gca,'YMinorTick','on');
    end
    drawnow;
    
    % stop case
    if all(abs(localMax_I-A_desired) <= intensityError) && ...
            all(abs(localMax_Ez(eff_ind)-Ez_desired(eff_ind)) <= intensityError)
        disp('Intensity is OK! Time Reversal Finished!');
        break;
    end
    
    if iteStep > 20
        if (std(STD_I(end-20:end)) < stdThreshold)
            disp('Standard Deviation is OK! Time Reversal Finished!');
            break;
        end
    end
end
%% show pupil
figure;
subplot(131);
pupilshowplr(plr,11);
subplot(132);
pupilshow(amp);colorbar;
subplot(133);
pupilshow(phs);colorbar;