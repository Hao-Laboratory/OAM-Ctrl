function out = MatrixTransfer(theta,phi,in,theta_max)

[M,N] = size(in);
xx = -1:2/(N-1):1;
yy = -1:2/(M-1):1;

% Becuase ( rho = h ), h = f*sin(theta),so h_step / h_max = sin(theta_step) / sin(theta_max)
% that is to say, rho_step / rho_max = sin(theta_step) / sin(theta_max) = [0:1]
xi = sin(theta)/sin(theta_max).*cos(phi);
yi = sin(theta)/sin(theta_max).*sin(phi);

out = interp2(xx,yy,in,xi,yi,'linear');