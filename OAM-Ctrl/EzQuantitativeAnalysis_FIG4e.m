clear;
close all;

pupilDiaPixNum = 300;
pixelSize = 0.01;
ratioE2M = 0:0.01:1;
l = 1:1:6;

obj.NA = 0.95; obj.n = 1;
beam.wavelength = 1;
beam.abr = zeros(pupilDiaPixNum);
scope.xs = -2.5:pixelSize:2.5;
scope.ys = 0;
scope.zs = 0;

%% set dipoles
coord.X = 0;
coord.Y = 0;
coord.Z = 0;

Dipole.Coord.X = [coord.X, coord.X];
Dipole.Coord.Y = [coord.Y, coord.Y];
Dipole.Coord.Z = [coord.Z, coord.Z];

N = length(coord.X);
Dipole.Type = [repmat('e',1,N),repmat('m',1,N),];

Dipole.Psi.X = mat2cell(zeros(1,2*N),1,ones(1,2*N));
Dipole.Psi.Y = mat2cell(zeros(1,2*N),1,ones(1,2*N));

Dipole.A.X = zeros(1,2*N);
Dipole.A.Y = zeros(1,2*N);

%% start loop
IzMax = zeros(1,length(ratioE2M));
for jj = 1:length(l)
    phi = vortexPhasePlate(pupilDiaPixNum,l(jj));
        Dipole.Psi.Z = {phi,phi};
    for ii = 1:length(ratioE2M)
        Dipole.A.Z = [sqrt(ratioE2M(ii)),sqrt(1-ratioE2M(ii))];
        
        % time reversal
        [amp,phs,plr] = TimeReversal_vec(pupilDiaPixNum,beam.wavelength,obj.NA,obj.n,Dipole);
        beam.amp = amp;
        beam.phs = phs;
        beam.plr = plr;
        
        [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,pupilDiaPixNum,...
            1,500,1);
        Ix = abs(Ex).^2; Iy = abs(Ey).^2; Iz = abs(Ez).^2; PSF = Ix + Iy + Iz;
        
        % normalization
        PSF_norm = PSF./max(PSF(:));
        Iz_norm = Iz./max(PSF(:));
        
        IzMax(ii) = max(Iz_norm(:));
    end
figure(1);
plot( ratioE2M,IzMax,'DisplayName',['\it l = ',num2str(l(jj))] );
axis tight;grid on;hold on;
ylim([0 0.8]);
legend;
set(gca,'box','off','XMinorTick','on','YMinorTick','on');
end
