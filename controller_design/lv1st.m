function [A,B] = lv1st(curr_x,curr_u,param)
% % x represents states as follows nx1
% % only designed for one single state
% 1:3 - r, position of the launch vehicle, iner
% 4:6 - v, velocity of the launch vehicle, iner
% 7:10 - q, quaternion from iner frame to body frame
% 11:13 - omega, angular velocity of the body, iner
% 14 - m, mass of the launch vehicle
% % % the following states are neglected for now
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
% 3:6 (12:15) - Grid fin actuation angle or angle of attack
% % % the following controls are neglected for now
% % 3 - "throttle" of main engine
% % 4:11 - on-off of the ACS thrusters

% decode dynamic states
r = curr_x(1:3,:);
v = curr_x(4:6,:);
q = curr_x(7:10,:);
omega = curr_x(11:13,:);
m = curr_x(14,:);
mox = curr_x(15,:);
mfuel = curr_x(16,:);

% decode controls
tvc_angle = curr_u(1:2,:);
gf_aoa = curr_u(3:6,:);
y = tvc_angle(1); p = tvc_angle(2);

% parameters
rn = vecnorm(r);
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
wx = omega(1); wy = omega(2); wz = omega(3);
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

% moment arms
rcp = param.rcp_1st;
renge = [-param.d.enge1st; 0; 0];
rc = param.rc1st(m,mox,mfuel);
moic = param.moic1st(m,mox,mfuel);
moicinv = inv(moic);
rcg2cp = rcp - rc;
rcg2eg = renge - rc;

% % % Compute state matrix A % % %
% Jacobian - rdot terms
Arr = zeros(3);
Arv = eye(3);
Arq = zeros(3,4);
Arw = zeros(3);
Arm = zeros(3,1);

% Jacobian - vdot terms
Avr_grav = 3*mu/((r'*r)^2.5)*r*r' - mu/((r'*r)^1.5)*eye(3);
Avr_aero = b/(2*m)*(vreln*OMEGA - e_d*vrel'*OMEGA);
dotCmat = [q0 q3 -q2; q1 q2 q3; -q2 q1 -q0; -q3 q0 q1;
    -q3 q0 q1; q2 -q1 q0; q1 q2 q3; -q0 -q3 q2;
    q2 -q1 q0; q3 -q0 -q1; q0 q3 -q2; q1 q2 q3]';

Avr = Avr_grav + Avr_aero;
Avv = -b/(2*m)*(vreln*eye(3) - e_d*vrel');
Avq = 2/m*Tmax*reshape(tvc_body'*dotCmat,[4,3])';
Avw = zeros(3);
Avm = -1/m^2*(Ft+Fd);

% Jacobian - qdot terms
Aqr = zeros(4,3);
Aqv = zeros(4,3);
Aqq = [0 -wx -wy -wz; wx 0 wz -wy; wy -wz 0 wx; wz wy -wx 0];
Aqw = [-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0];
Aqm = zeros(4,1);

% Jacobian - wdot terms
dotMmatr = [cross(rcg2cp,Avr_aero(:,1)) cross(rcg2cp,Avr_aero(:,2)) cross(rcg2cp,Avr_aero(:,3))]*m;
dotMmatv = [cross(rcg2cp,Avv(:,1)) cross(rcg2cp,Avv(:,2)) cross(rcg2cp,Avv(:,3))]*m;
dotMmatq = [cross(rcg2eg,Avq(:,1)) cross(rcg2eg,Avq(:,2)) cross(rcg2eg,Avq(:,3)) cross(rcg2eg,Avq(:,4))]*m;
dotWmatw = [-moic(3,1)*wy+moic(2,1)*wz moic(3,2)*wx-moic(1,2)*wz -moic(2,3)*wx+moic(1,3)*wy;
    -2*moic(2,1)*wx+(moic(1,1)-moic(2,2))*wy-moic(2,3)*wz (moic(1,1)-moic(2,2))*wx+2*moic(1,2)*wy+moic(1,3)*wz (moic(3,3)-moic(1,1))*wx-moic(1,2)*wy-2*moic(1,3)*wz;
    2*moic(3,1)*wx+moic(3,2)*wy+(moic(3,3)-moic(1,1))*wz -moic(3,1)*wx-2*moic(3,2)*wy+(moic(2,2)-moic(3,3))*wz moic(2,1)*wx+(moic(2,2)-moic(3,3))*wy+2*moic(2,3)*wz];

Awr = moicinv*dotMmatr;
Awv = moicinv*dotMmatv;
Awq = moicinv*dotMmatq;
Aww = moicinv*dotWmatw;
Awm = zeros(3,1);

% Jacobian - mdot terms
Ftn = vecnorm(Ft);
cstar = param.g0*param.Isp1;

Amr = zeros(1,3);
Amv = zeros(1,3);
Amq = Ft'*Avq*m/(cstar*Ftn);
Amw = zeros(1,3);
Amm = 0;

A = [Arr Arv Arq Arw Arm;
    Avr Avv Avq Avw Avm;
    Aqr Aqv Aqq Aqw Aqm;
    Awr Awv Awq Aww Awm;
    Amr Amv Amq Amw Amm];

% % % Compute input matrix B % % %
% Jacobian - rdot terms
Bry = zeros(3,1);
Brp = zeros(3,1);

% Jacobian - vdot terms
Bvy = C*Tmax*[-cos(p)*sin(y); cos(p)*cos(y); 0]/m;
Bvp = C*Tmax*[-sin(p)*cos(y); -sin(p)*sin(y); cos(p)]/m;

% Jacobian - qdot terms
Bqy = zeros(4,1);
Bqp = zeros(4,1);

% Jacobian - wdot terms
Bwy = moicinv*cross(rcg2eg,Bvy)*m;
Bwp = moicinv*cross(rcg2eg,Bvp)*m;

% Jacobian - mdot terms
Bmy = Ft'*Bvy*m/(cstar*Ftn);
Bmp = Ft'*Bvp*m/(cstar*Ftn);

B = [Bry Brp;
    Bvy Bvp;
    Bqy Bqp;
    Bwy Bwp;
    Bmy Bmp];

end