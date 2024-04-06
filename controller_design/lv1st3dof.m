function [A,B] = lv1st3dof(x,u,param)
% % x represents states as follows nx1
% % only designed for one single state
% 1:3 - r, position of the launch vehicle, iner
% 4:6 - v, velocity of the launch vehicle, iner
% % 7:10 - q, quaternion from iner frame to body frame
% % 11:13 - omega, angular velocity of the body, iner
% 7(14) - m, mass of the launch vehicle
% % 15 - mox, mass of the oxidizer
% % 16 - mfuel, mass of the fuel
% % 17:18 - phiox and thetaox, body
% %   phiox, oxidizer sloshing mass off-b2 pitch anlge
% %   thetaox, oxidizer sloshing mass off-b3 roll angle
% % 19:20 - phifuel and thetafuel, body
% %   phifuel, fuel sloshing mass off-b2 pitch anlge
% %   thetafuel, fuel sloshing mass off-b3 roll angle
% % 21:22 - dphiox and dthetaox, body
% % 23:24 - dphifuel and dthetafuel, body

% % u represents control as follows mx1
% 1:2 - TVC yaw and pitch angles, body
% % % the following controls are neglected for now
% % 3 - "throttle" of main engine
% % 4:11 - on-off of the ACS thrusters

% decode dynamic states
r = x(1:3,:);
v = x(4:6,:);
m = x(14,:);
q = x(7:10,:);

% decode controls
tvc_angle = u(1:2,:);
y = tvc_angle(1); p = tvc_angle(2);

% parameters
rn = vecnorm(r);
mu = param.mu;
Tmax = param.maxT_1st;

% regulate quaternion
q = q./vecnorm(q);
C = EP2C(q);

% TVC
tvc_body = param.tvc_point(y,p);
tvc_iner = C*tvc_body;
Ft = Tmax*tvc_iner;

% atmosphere
h = rn-param.earthR;
atmo_profile = atmo(h);
rho = atmo_profile.rho;
Cd = 0.5*ones(size(rho));
b = Cd.*rho.*param.S;

% aerodynamic drag
OMEGA = param.skewOMEGA;
atmo_v = OMEGA*r;
vrel = v - atmo_v;
vreln = vecnorm(vrel);
e_d = -vrel./vreln;
Fd = 1/2.*b.*vreln.^2.*e_d;

% % % Compute state matrix A % % %
% Jacobian - rdot terms
Arr = zeros(3);
Arv = eye(3);
Arm = zeros(3,1);

% Jacobian - vdot terms
Avr_grav = 3*mu/((r'*r)^2.5)*r*r' - mu/((r'*r)^1.5)*eye(3);
Avr_aero = b/(2*m)*(vreln*OMEGA - e_d*vrel'*OMEGA);

Avr = Avr_grav + Avr_aero;
Avv = -b/(2*m)*(vreln*eye(3) - e_d*vrel');
Avm = -1/m^2*(Ft+Fd);

% Jacobian - mdot terms
Ftn = vecnorm(Ft);
cstar = param.g0*param.Isp1;

Amr = zeros(1,3);
Amv = zeros(1,3);
Amm = 0;

A = [Arr Arv Arm;
    Avr Avv Avm;
    Amr Amv Amm];

% % % Compute input matrix B % % %
% Jacobian - rdot terms
Bry = zeros(3,1);
Brp = zeros(3,1);

% Jacobian - vdot terms
Bvy = C*Tmax*[-cos(p)*sin(y); cos(p)*cos(y); 0]/m;
Bvp = C*Tmax*[-sin(p)*cos(y); -sin(p)*sin(y); cos(p)]/m;

% Jacobian - mdot terms
Bmy = Ft'*Bvy*m/(cstar*Ftn);
Bmp = Ft'*Bvp*m/(cstar*Ftn);

B = [Bry Brp;
    Bvy Bvp;
    Bmy Bmp];

end