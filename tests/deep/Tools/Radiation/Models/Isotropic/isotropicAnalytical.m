% Verification of the analytical results for the
% isotropic model in SOFTv2.
clc; clear;

x = [0.75,0.2,0.04]';       % Particle position
X0 = [0,-1.069,0]';         % Detector position
l = 0.006;                  % Detector aperture

xhat = [1,0,0]';
yhat = [0,1,0]';
zhat = [0,0,1]';

nhat = yhat;
ehat1 = zhat;
ehat2 = xhat;

cn = nhat' * (x - X0);
c1 = ehat1' * (x - X0);
c2 = ehat2' * (x - X0);

intg = @(lx,ly) 1./(cn^2 + (c1-lx).^2 + (c2-ly).^2).^1.5;

In = cn*integral2(@(x,y) intg(x,y), -l, l, -l, l);

% Analytical result
Liso = @(x,y) atan(((c1-x).*(c2-y)) ./ (cn.*sqrt(cn^2 + (c2-y).^2 + (c1-x).^2)));

Ia = Liso(l,l) + Liso(-l,-l) - Liso(-l,l) - Liso(l,-l);

disp(['Numerical:   ',num2str(In)]);
disp(['Analytical:  ',num2str(Ia)]);
disp(['Delta:       ',num2str(abs((In-Ia)/In)*100),'%']);