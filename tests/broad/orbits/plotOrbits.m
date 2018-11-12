clc; clear;

load gc-orbit
%load gc-orbit-drifts

nt = numel(t);
SOL = reshape(solution, [6,nt]);
X = reshape(x, [3,nt]);
P = reshape(p, [3,nt]);
XX = SOL(1:3,:);
PP = SOL(4:6,:);

R = sqrt(X(1,:).^2 + X(2,:).^2);
Z = X(3,:);
RR = sqrt(XX(1,:).^2 + XX(2,:).^2);
ZZ = XX(3,:);

figure(1);
set(gcf, 'Position', [500, 400, 400, 250]);
box on
hold on
plot(R, Z, 'linewidth', 2);
plot(RR, ZZ);
xlabel('$R$ (m)', 'Interpreter', 'latex');
ylabel('$Z$ (m)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
title('Orbit');
axis equal

figure(2);
set(gcf, 'Position', [1000, 400, 400, 250]);
box on
hold on
PMAG = sqrt(P(1,:).^2 + P(2,:).^2 + P(3,:).^2);
PPMAG = sqrt(PP(1,:).^2 + PP(2,:).^2 + PP(3,:).^2);
plot(t, PMAG);
%plot(t, PPMAG, '--');
xlabel('$t$ (s)', 'Interpreter', 'latex');
ylabel('$p$ ($mc$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
title('Momentum');
