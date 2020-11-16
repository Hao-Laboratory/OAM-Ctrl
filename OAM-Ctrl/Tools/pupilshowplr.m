function pupilshowplr(plr,num,varargin)
%PUPILSHOWPLR show the polarization distribution on the pupil plane
%   
% author: Xin Liu
% email:  liuxin2018@zju.edu.cn
% Oct.14, 2020

switch nargin
    case 2
        Color = 'k';
    case 3
        Color = varargin{1};
end

px = plr(:,:,1);
py = plr(:,:,2);
% pz = plr(:,:,3);  % pz == 0.

scale = round(length(px)/num);
px = px(1:scale:length(px),1:scale:length(px));
py = py(1:scale:length(py),1:scale:length(py));
N = length(px);

Ax1 = abs(px);
Ay1 = abs(py);
Ax = Ax1./sqrt(Ax1.^2+Ay1.^2);
Ay = Ay1./sqrt(Ax1.^2+Ay1.^2);
phs_x = angle(px);
phs_y = angle(py);

delta_phs = phs_y-phs_x;

[xx, yy] = meshgrid(linspace(-N+1,N-1,N),linspace(-N+1,N-1,N));
[~,rho] = cart2pol(xx,yy);
Ax(rho>N-1) = nan;
Ay(rho>N-1) = nan;

phi = linspace(0,2*pi,361);

% arrows parameters
h = 0.2; % height
w = 0.2; % width
ar = [-w/2,h; 0,0; w/2,h];

r = 0.66;

case1 = abs(sin(delta_phs)) < 1e-5;  % delta phase is m*pi
case2 = abs(cos(delta_phs)) < 1e-5;  % delta phase is (2m+1)*pi/2
case3 = mod(delta_phs,2*pi) < pi;  % y delay x

%% k1
k1 = find(case1);  % linear polarization

quiver(xx(k1),yy(k1),sign(cos(phs_x(k1))).*Ax(k1), sign(cos(phs_y(k1))).*Ay(k1),...
    Color,'MaxHeadSize',1,'AutoScaleFactor',0.5);
hold on;

%% k
% elliptical or circular polarization
k2 = find(case2 & case3);  % right hand circular
k3 = find(case2 & ~case3);  % left hand circular
k4 = find(~case1 & ~case2 & case3);  % right hand elliptical
k5 = find(~case1 & ~case2 & ~case3);  % left hand elliptical

alpha_k2 = 0;
alpha_k3 = 0;

alpha_k4 = 0.5*atan( 2*Ax(k4).*Ay(k4).*cos(delta_phs(k4)) ./ (Ax(k4).^2-Ay(k4).^2) );
alpha_k5 = 0.5*atan( 2*Ax(k5).*Ay(k5).*cos(delta_phs(k5)) ./ (Ax(k5).^2-Ay(k5).^2) );

arx_k2 = repmat(ar(:,1),1,numel(k2));
ary_k2 = repmat(-ar(:,2),1,numel(k2));

arx_k3 = repmat(ar(:,1),1,numel(k3));
ary_k3 = repmat(ar(:,2),1,numel(k3));

arx_k4 = repmat(ar(:,1),1,numel(k4));
ary_k4 = repmat(-ar(:,2),1,numel(k4));

arx_k5 = repmat(ar(:,1),1,numel(k5));
ary_k5 = repmat(ar(:,2),1,numel(k5));

ax_k2 = arx_k2+r.*Ax(k2).';
ay_k2 = ary_k2;

ax_k3 = arx_k3+r.*Ax(k3).';
ay_k3 = ary_k3;

ax_k4 = arx_k4+r.*Ax(k4).';
ay_k4 = ary_k4;

ax_k5 = arx_k5+r.*Ax(k5).';
ay_k5 = ary_k5;

%% k2 
% the ellipse at the origin
xo_k2 = r*Ax(k2).*cos(phi);
yo_k2 = r*Ay(k2).*sin(phi);
xo_k3 = r*Ax(k3).*cos(phi);
yo_k3 = r*Ay(k3).*sin(phi);
xo_k4 = r*Ax(k4).*cos(phi);
yo_k4 = r*Ay(k4).*sin(phi);
xo_k5 = r*Ax(k5).*cos(phi);
yo_k5 = r*Ay(k5).*sin(phi);

% plot the ellipse
x2_k2 = xx(k2) + xo_k2.*cos(alpha_k2) - yo_k2.*sin(alpha_k2);
y2_k2 = yy(k2) + xo_k2.*sin(alpha_k2) + yo_k2.*cos(alpha_k2);
x2_k3 = xx(k3) + xo_k3.*cos(alpha_k3) - yo_k3.*sin(alpha_k3);
y2_k3 = yy(k3) + xo_k3.*sin(alpha_k3) + yo_k3.*cos(alpha_k3);
x2_k4 = xx(k4) + xo_k4.*cos(alpha_k4) - yo_k4.*sin(alpha_k4);
y2_k4 = yy(k4) + xo_k4.*sin(alpha_k4) + yo_k4.*cos(alpha_k4);
x2_k5 = xx(k5) + xo_k5.*cos(alpha_k5) - yo_k5.*sin(alpha_k5);
y2_k5 = yy(k5) + xo_k5.*sin(alpha_k5) + yo_k5.*cos(alpha_k5);

plot(x2_k2.', y2_k2.', Color);hold on;
plot(x2_k3.', y2_k3.', Color);hold on;
plot(x2_k4.', y2_k4.', Color);hold on;
plot(x2_k5.', y2_k5.', Color);hold on;

% plot the arrow
arrowX_k2 = ax_k2.*cos(alpha_k2) - ay_k2.*sin(alpha_k2);
arrowY_k2 = ax_k2.*sin(alpha_k2) + ay_k2.*cos(alpha_k2);

arrowX_k3 = ax_k3.*cos(alpha_k3) - ay_k3.*sin(alpha_k3);
arrowY_k3 = ax_k3.*sin(alpha_k3) + ay_k3.*cos(alpha_k3);

arrowX_k4 = ax_k4.*cos(alpha_k4.') - ay_k4.*sin(alpha_k4.');
arrowY_k4 = ax_k4.*sin(alpha_k4.') + ay_k4.*cos(alpha_k4.');

arrowX_k5 = ax_k5.*cos(alpha_k5.') - ay_k5.*sin(alpha_k5.');
arrowY_k5 = ax_k5.*sin(alpha_k5.') + ay_k5.*cos(alpha_k5.');

arrowX_k2 = xx(k2).'+arrowX_k2;
arrowY_k2 = yy(k2).'+arrowY_k2;
arrowX_k3 = xx(k3).'+arrowX_k3;
arrowY_k3 = yy(k3).'+arrowY_k3;
arrowX_k4 = xx(k4).'+arrowX_k4;
arrowY_k4 = yy(k4).'+arrowY_k4;
arrowX_k5 = xx(k5).'+arrowX_k5;
arrowY_k5 = yy(k5).'+arrowY_k5;

plot(arrowX_k2,arrowY_k2,'r');hold on;
plot(arrowX_k3,arrowY_k3,'r');hold on;
plot(arrowX_k4,arrowY_k4,'r');hold on;
plot(arrowX_k5,arrowY_k5,'r');hold on;

hold off;axis tight;axis equal;axis off;
end
