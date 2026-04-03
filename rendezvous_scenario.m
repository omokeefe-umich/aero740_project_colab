% Rendezvous Scenario Simulation
clear; clc; close all;

% Construct a scenario where the quadcopter starts at a given state and 
% must rendezvous with a target vehicle before it achieves a certain position.

% The target vehicle is moving along a known trajectory, and the quadcopter must 
% adjust its control inputs to meet the target prior to a specific point in time.

% The vehicle is assumed to at an altitude of 0m and moving at a speed of Vt = 5 m/s 
% along the x-axis. The quadcopter starts at an altitude of 10m and must descend to 
% meet the target before it is capable of reaching a position of x = 50m. The 
% rendezvous point is defined as the point where the quadcopter's position 
% matches the target's position within a certain tolerance.

% Vehicle State vectors (12x1):
%   [x, y, z, phi, theta, psi, xdot, ydot, zdot, phidot, thetadot, psidot]
%   positions (m), Euler angles (rad), linear velocities (m/s), angular rates (rad/s)

% Ground Vehicle Trajectory
Vmax = 10; % m/s
d0 = 100; % Starting positon in the reference frame
d_gate = 350; % The x-distance of the gate in the reference frame
t_gate = d_gate / Vmax; % Time at which the target reaches the gate
d_target = @(t) d0 + Vmax * t; % Target's position as a function of time

x_gv = [d_target(0); 0; 0; 0; 0; 0; Vmax; 0; 0; 0; 0; 0]; % Ground vehicle state vector

% Quadcopter
m_kg  = 1.12;        % mass [kg]
g_ms2  = 9.81;        % gravitational acceleration [m/s^2]
l_m  = 0.23;        % arm length [m]
Jr_kgm2 = 8.5e-4;      % rotor moment of inertia [kg.m^2]
Ix_kgm2 = 0.0019;      % roll moment of inertia [kg.m^2]
Iy_kgm2 = 0.0019;      % pitch moment of inertia [kg.m^2]
Iz_kgm2 = 0.0223;      % yaw moment of inertia [kg.m^2]
b_nrads2  = 7.73212e-6;  % thrust coefficient [N/(rad/s)^2]
d_nrads2  = 1.27513e-7;  % drag (torque) coefficient [N.m/(rad/s)^2]

% Initial state of the quadcopter (12x1)
x0 = [0; 0; 100; 10; 0; -5; 0; 0; 0; 0; 0; 0]; 
u0 = [m*g; 0; 0; 0]; % Initial control input vector (4x1)

% time discretization
N = 10;
t_final = t_gate; % Final time for the simulation
dt = t_final / N; % Time step

% Discretized Quadcopter Dynamics (using a simple Euler integration for demonstration)
A_d = [ 1.0000, 	0, 	0, 	0    0.0491, 	0    0.1000, 	0, 	0, 	0    0.0016, 	0;
        0    1.0000, 	0   -0.0491, 	0, 	0, 	0    0.1000, 	0   -0.0016, 	0, 	0;
        0, 	0    1.0000, 	0, 	0, 	0, 	0, 	0    0.1000, 	0, 	0, 	0;
        0, 	0, 	0    1.0000, 	0, 	0, 	0, 	0, 	0    0.1000, 	0, 	0;
        0, 	0, 	0, 	0    1.0000, 	0, 	0, 	0, 	0, 	0    0.1000, 	0;
        0, 	0, 	0, 	0, 	0    1.0000, 	0, 	0, 	0, 	0, 	0    0.1000;
        0, 	0, 	0, 	0    0.9810, 	0    1.0000, 	0, 	0, 	0    0.0491, 	0;
        0, 	0, 	0   -0.9810, 	0, 	0, 	0    1.0000, 	0   -0.0491, 	0, 	0;
        0, 	0, 	0, 	0, 	0, 	0, 	0, 	0    1.0000, 	0, 	0, 	0;
        0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0    1.0000, 	0, 	0;
        0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0    1.0000, 	0;
        0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	0    1.0000];

B_d = [0, 	0    0.0016, 	0;
 	0   -0.0020, 	0, 	0;
   0.0033, 	0, 	0, 	0;
 	0    0.2500, 	0, 	0;
 	0, 	0    0.2000, 	0;
 	0, 	0, 	0    0.1667;
 	0, 	0    0.0654, 	0;
 	0   -0.0818, 	0, 	0;
   0.0667, 	0, 	0, 	0;
 	0    5.0000, 	0, 	0;
 	0, 	0    4.0000, 	0;
 	0, 	0, 	0    3.3333];

x_qc = x0; % Initialize quadcopter state
u = u0; % Initialize control input
x_qc_kp1 = A_d * x_qc + B_d * u;

function dx = combined_dynamics(t, x, A_d, B_d, u, x_gv)
    dx = A_d * x + B_d * u - x_gv; % Simple combined dynamics for demonstration
end