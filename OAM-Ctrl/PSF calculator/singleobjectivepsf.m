function [Ex,Ey,Ez] = singleobjectivepsf(obj,beam,scope,...
    pupilDiaPixNum)

% SINGLEOBJECTIVEPSF calculates the PSF of normal microscopy with one
% single objective.
%
% To calculate the result. we need to define the amplitude (amp),
% polarization (px, py, pz), aberration (abr) and phase (phs) are identical
% for the incident beam. We further define the numerical aperture (NA) and
% the refractive index (n) for the objective.
%
% The calculation scale is defined by xs, ys, and zs along different axes
% in 3D.
%
% format: [Ex,Ey,Ez] = SingleObjectivePSF(NA,n,abr,amp,px,py,pz,phs,...
%                      xs,ys,zs)
%
% Note for input:
% NA:         single value
% n:          single value, n > NA
% abr:        2D square matrix with the same size (m*m)
% amp:        2D square matrix with the same size (m*m)
% px,py,pz:   2D square matrix with the same size (m*m)
% phs:        2D square matrix with the same size (m*m)
% xs:         1D array
% ys:         1D array
% zs:         1D array
%
% coded by HAO,Xiang
% xiang.hao@yale.edu
% Jan 26, 2017
%
% final updated by Xin Liu
% Feb, 2019

NA = obj.NA;
n = obj.n;
lambda = beam.wavelength;
abr = beam.abr;
amp = beam.amp;
px = beam.plr(:,:,1);
py = beam.plr(:,:,2);
pz = beam.plr(:,:,3);
phs = beam.phs;
xs = scope.xs;
ys = scope.ys;
zs = scope.zs;

%% input data preparation
k = 2*pi/lambda; % wave vector
theta_max = asin(NA/n); % maximum divergent angle
theta_number = pupilDiaPixNum; % determine the calculation speed and accuracy
phi_number = pupilDiaPixNum; % determine the calculation speed and accuracy
theta_step = 0:asin(NA/n)/theta_number:asin(NA/n);
phi_step = 0 : 2*pi/phi_number : 2*pi;
[theta,phi] = meshgrid(theta_step,phi_step);
ampT = MatrixTransfer(theta,phi,amp,theta_max);
pxT = MatrixTransfer(theta,phi,px,theta_max);
pyT = MatrixTransfer(theta,phi,py,theta_max);
pzT = MatrixTransfer(theta,phi,pz,theta_max);
plrT = cat(3,pxT,pyT,pzT);
phsT = MatrixTransfer(theta,phi,phs,theta_max);
abrT = MatrixTransfer(theta,phi,abr,theta_max);

%% calculation initialization
% u and l correspond to upper and lower beams, respectively
lx = length(xs);  ly = length(ys);  lz = length(zs);
flag = [(lx == 1),(ly == 1),(lz == 1)];     % calculation flag
flag = bi2de(flag,'left-msb');
Ex = zeros(ly,lx,lz);
Ey = Ex;
Ez = Ex;

%% calculation
totalComp = lx * ly * lz;   % total number of components to be calculated

dispstat('','init');
for j = 1:lz
    z2 = zs(j);
    for jj = 1:ly
        for jjj = 1:lx
            [phi2,r2] = cart2pol(xs(jjj),ys(jj));
            
            [dEx,dEy,dEz] = debye2(theta,phi,r2,z2,phi2,abrT,ampT,plrT,...
                phsT,k,n);
            
            % calculate Ex,Ey,Ez
            Ex(jj,jjj,j) = trapez2(theta_step,phi_step,dEx);
            Ey(jj,jjj,j) = trapez2(theta_step,phi_step,dEy);
            Ez(jj,jjj,j) = trapez2(theta_step,phi_step,dEz);
            
            counter = jjj + (jj - 1) * lx + (j - 1) * ly * lx;
            dispstat([num2str(counter/totalComp*100),' %']);
        end
    end
end
dispstat('','keepprev');

%% transpose the matrix
switch flag
    case 0   % xyz mode
        disp('XYZ calculation complete!');
    case 1   % xy mode
        Ex = squeeze(Ex);
        Ey = squeeze(Ey);
        Ez = squeeze(Ez);
        disp('XY calculation complete!');
    case 2   % xz mode
        Ex = squeeze(Ex);    Ex = rot90(Ex);   Ex = flipud(Ex);
        Ey = squeeze(Ey);    Ey = rot90(Ey);   Ey = flipud(Ey);
        Ez = squeeze(Ez);    Ez = rot90(Ez);   Ez = flipud(Ez);
        disp('XZ calculation complete!');
    case 3   % x mode
        disp('X calculation complete!');
    case 4   % yz mode
        Ex = squeeze(Ex);    Ex = rot90(Ex);   Ex = flipud(Ex);
        Ey = squeeze(Ey);    Ey = rot90(Ey);   Ey = flipud(Ey);
        Ez = squeeze(Ez);    Ez = rot90(Ez);   Ez = flipud(Ez);
        disp('YZ calculation complete!');
    case 5   % y mode
        Ex = squeeze(Ex);
        Ey = squeeze(Ey);
        Ez = squeeze(Ez);
        disp('Y calculation complete!');
    case 6   % z mode
        Ex = squeeze(Ex);
        Ey = squeeze(Ey);
        Ez = squeeze(Ez);
        disp('Z calculation complete!');
    case 7   % point mode
        warning('point calculation complete!');
    otherwise
        error('calculation flag invalid! check input!');
end