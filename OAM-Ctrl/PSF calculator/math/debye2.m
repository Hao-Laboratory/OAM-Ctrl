function [dx,dy,dz]=debye2(theta,phi,r2,z2,phi2,abr,amp,plr,phs,k,n)
% DEBYE2 is the Debye integral expression but using matrix as input.
%
% Corresponding expression is fundamentally from B. Richards and E. Wolf, 
% "Electromagnetic Diffraction in Optical Systems .2. Structure of the 
% Image Field in an Aplanatic System," Proc R Soc Lon Ser-A 253, 358-379 
% (1959).
%
% input:
% *************************************************
% theta: flare angle
% phi: azimuthal angle
% r2,z2,phi2: polar coordinate that defines the position in the focal field
% abe: aberration matrix
% amp: amplitude matrix
% plr: polarization matrix
% phs: phase deley matrix
% k: beam vector
% n: refractive index
%
% coded by HAO,Xiang
% first coded on Jul. 25, 2014
% last updated on Jul. 25, 2014

prefix = sin(theta);
suffix = exp(1i*k*n*(z2*cos(theta)-r2*sin(theta).*cos(phi-phi2)));
lens = AL(theta);  % aplantic lens
aberration = exp(1i*abr);  % aberration
phase = exp(1i*phs);  % phase
% polarization status corresponding to Ex, Ey and Ez
Px = (1+(cos(theta)-1).*((cos(phi)).^2)).*plr(:,:,1)+ ...
    ((cos(theta)-1).*cos(phi).*sin(phi)).*plr(:,:,2)+ ...
    -sin(theta).*cos(phi).*plr(:,:,3);
Py = (cos(theta)-1).*cos(phi).*sin(phi).*plr(:,:,1)+ ...
    (1+(cos(theta)-1).*(sin(phi)).^2).*plr(:,:,2)+ ...
    -sin(theta).*sin(phi).*plr(:,:,3);
Pz = sin(theta).*cos(phi).*plr(:,:,1)+ ...
    sin(theta).*sin(phi).*plr(:,:,2)+ ...
    cos(theta).*plr(:,:,3);

dx = prefix.*lens.*aberration.*amp.*Px.*phase.*suffix;
dy = prefix.*lens.*aberration.*amp.*Py.*phase.*suffix;
dz = prefix.*lens.*aberration.*amp.*Pz.*phase.*suffix;