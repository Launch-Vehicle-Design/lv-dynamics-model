clear; clc; close all
addpath("..\parameter_functions")

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

% environmental parameter
param.mu = 3.986004418e14;
param.earthR = 6357e3;
param.omega_len = 2*pi/(24*3600);
param.g = @(h) g(h);
param.rho = @(h) rho(h);
param.P = @(h) P(h);
param.gamma = 1.4;

% coordinate parameter
param.x = [1;0;0];
param.y = [0;1;0];
param.z = [0;0;1];

% launch vehicle parameter
param.isp_1st = 276; param.isp_2nd = 369.5;
% vehicle masses (subject to change)
param.ms_1st = 125.2035; param.ms_2nd = 29.1843;
param.mp_1st = 1801.0036; param.mp_2nd = 248.7618;
param.mpl_1st = 314.9462; param.mpl_2nd = 37;
% vehicle maximum thrust bounds
param.mtot_1st = param.ms_1st+param.mp_1st+param.mpl_1st;
param.mtot_2nd = param.ms_2nd+param.mp_2nd+param.mpl_2nd;
param.TtoW_1st = 1.6; param.TtoW_2nd = 0.7;
param.maxT_1st = param.TtoW_1st*param.mtot_1st*9.80665;
param.maxT_2nd = param.TtoW_2nd*param.mtot_2nd*9.80665;
% vehicle length assuming no COM movement (subject to change)
% -------------------------------------------------------------
% 1st stage solid | interstage | 2nd stage RDRE | GNC | Payload >
% -------------------------------------------------------------
param.plf2eng_1st0 = -18*0.3048; param.plf2eng_2nd0 = -4.5*0.3048;
param.plf2cp_1st = -1.5*0.3048; param.plf2cp_2nd = -1.2*0.3048;

param.sitelat = 80; param.sitelon = -100;
% param.omegaE = omega_len*[0;sin(param.sitelat/180*pi);cos(param.sitelat/180*pi)];
param.vair = @(r) -param.x*param.omega_len*dot(param.y,r)*cos(param.sitelat);
param.orbinc = 100;
param.dia = 23.94*0.0254;
param.S = pi*param.dia^2/4;
param.CL = @(M,AoA) CL(M,AoA);
param.CD = @(M,AoA) CD(M,AoA);

function [dxdt,stage] = rocketDyn(t,x,u,stage,param)
% inertial reference frame based dynamics
%
% x - dynamic state that monitors the launch vehicle's state
% 1:3 - position in the sliced ECI frame by launch location
% 4:6 - velocity in the sliced ECI frame by launch location
% 7 - pitch angle
% 8 - angular velocity of the pitch angle
% 9 - vehicle mass
% 
% u - control input that drives launch vehicle's orbit insertion
% 1 - delta thrust angle (can be achieved by different TVCs)
% 2 - throttle percentage of nomimal maximum thrust
%
% stage - 0: lv drop; 1: 1st stg; 2: 2nd stg

% state transcripts
% r and v here should be pseduo-3D where rz and vz are 0s
r = x(1:3); v = x(4:6);
rlen = norm(r); vlen = norm(v);
iner_pa = x(7); pitch_dot = x(8);
m = x(9);
% control transcripts
delta = u(1); throttle = u(2);

% non-inertial coordinates
e_lv = r/rlen; e_lh = cross(r,param.z)/rlen;
iner2body = [cos(iner_pa) -sin(iner_pa) 0; sin(iner_pa) cos(iner_pa) 0; 0 0 1];
e_by = iner2body*param.y; e_bx = iner2body*param.x;

% derived parameters
h = rlen-param.earthR;
a = sqrt(param.gamma*param.P(h)/param.rho(h));
pa = acos(dot(e_by,e_lh))*sign(cross(e_lh,e_by));

% apply relative atmosphere flow in inertial frame
vair = param.vair(r); vrela = v - vair;
vrelalen = norm(vrela); mach = vrelalen/a;
aoa = acos(dot(vrela,e_by)/vlen)*sign(cross(e_by,vrela));
% define lift-drag coordinate frame
e_drag = vrela/vrelalen;
e_lift = cross(param.z,vrela)/norm(cross(param.z,vrela));

% external force profile in 2D inertial frame
Fg = -param.mu*m/(rlen^3)*r;
FD = -1/2*param.rho(h)*param.S*param.CD(mach,aoa*180/pi)*vrelalen^2*e_drag;
FL = 1/2*param.rho(h)*param.S*param.CL(mach,aoa*180/pi)*vrelalen^2*e_lift;
% stage dependant forces
Fthrust = 0; mdot = 0;
if stage == 1
    Fthrust = throttle*param.maxT_1st*[cos(pa+delta);sin(pa+delta);0];
    mdot = Fthrust/(param.g(0)*param.isp_1st);
elseif stage == 2
    Fthrust = throttle*param.maxT_2nd*[cos(pa+delta);sin(pa+delta);0];
    mdot = Fthrust/(param.g(0)*param.isp_2nd);
end
% moment resulted from external moment reference to payload fairing tip


dxdt = zeros(size(x));
% kinematics in inertial frame
dxdt(1:3) = v;
% Newton 2nd law in inertial frame
dxdt(4:6) = 1/m*(Fg + FD + FL + Fthrust);
% rotational kinematics in inertial frame
dxdt(7) = pitch_dot;
% Euler 2nd law in non-CG non-inertial frame
dxdt(8) = 0;
dxdt(9) = -mdot;
end