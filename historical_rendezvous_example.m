clear all; clc;close all;
P0 = [0 0];
E0 = [100 0];
VE = 5; 
VP = 6;
thetaE= pi/2;
R0 = norm(E0 - P0);
beta0 = atan2(E0(2) - P0(2), E0(1) - P0(1));
x0 = [R0; beta0];
tspan = [0:0.01:100];
opts = odeset('Events', @stopFcn, 'RelTol', 1e-6, 'AbsTol', 1e-6);
[t, x] = ode45(@(t,x) RandD(t,x,VP,VE,thetaE), tspan, x0, opts);
R = x(:,1);
beta = x(:,2);
xE = E0(1) + VE*cos(thetaE).*t;
yE = E0(2) + VE*sin(thetaE).*t;
xP = xE - R.*cos(beta);
yP = yE - R.*sin(beta);
Tc_th = R0 / (VP - VE);
Tc_num = t(end);

%% Plots
figure; hold on; axis equal; grid on; box on;
plot(xE, yE, 'r-', 'LineWidth', 1.5, 'DisplayName','Evader');
plot(xP, yP, 'b-', 'LineWidth', 1.5, 'DisplayName','Pursuer');
plot(P0(1), P0(2), 'bo', 'MarkerFaceColor','b', 'DisplayName','P0');
plot(E0(1), E0(2), 'ro', 'MarkerFaceColor','r', 'DisplayName','E0');
legend('Location','best');
xlabel('x [m]');
ylabel('y [m]');
% Range vs time
figure; plot(t,R,'LineWidth',1.5);
grid on; box on;
xlabel('t [s]');
ylabel('R(t) [m]');
title('Range vs Time');
%%
function dx = RandD(~, x, VP, VE, thetaE)
    R = x(1);
    beta = x(2);
    dR = VE*cos(beta - thetaE) - VP;
    db = -VE*sin(beta - thetaE)/R;
    dx = [dR; db];
end
function [value, isterminal, direction] = stopFcn(~, x)
    R = x(1);
    value = R;
    isterminal = 1;
    direction = -1;
end