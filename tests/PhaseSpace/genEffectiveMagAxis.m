% Generate a table of effective magnetic-axis values
%
% Written by: Mathias Hoppe, hoppe@chalmers.se
% April 2018
clc; clear;

Rm = 0.68;
B0 = 5;

% Below, we assume theta = phi = 0 & q = 1
phihat = [0 1 0];
thetahat = [0 0 1];
rhat = [-1 0 0];
Babs = @(r) B0*Rm./(Rm-r) .* sqrt((r./Rm).^2+1);
B = @(r) B0*Rm ./ (Rm-r) .* (r/Rm * thetahat + phihat);
gradB = @(r) rhat * B0*Rm./(Rm-r) .* (r./(Rm^2*sqrt((r./Rm).^2+1)) +...
    1./(Rm-r).*sqrt((r/Rm).^2+1));
gradParBhat = @(r) (Rm-r)./(Rm+r) * rhat;
bhat = @(r) B(r) ./ Babs(r);

XdotPerp = @(r,ppar,pperp) 1./(g*B
