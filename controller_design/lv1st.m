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
y = tvc_angle(1);
p = tvc_angle(2);
cy = cos(y); sy = sin(y);
cp = cos(p); sp = sin(p);

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
cstar = param.g0*param.Isp1;
tvc_body = param.tvc_point(y,p);
tvc_iner = C*tvc_body;
Ftb = Tmax*tvc_body;
Ft = Tmax*tvc_iner;

% atmosphere
h = rn-param.earthR;
atmo_profile = atmo(h);
rho = atmo_profile.rho;
pres = atmo_profile.P;
a = sqrt(param.gamma.*pres./rho);
Cd = 0.5*ones(size(rho));
b = Cd.*rho.*param.S;

% aerodynamic drag
OMEGA = param.skewOMEGA;
atmo_v = OMEGA*r;
vrel = v - atmo_v;
vreln = vecnorm(vrel);
mach = vreln./a;
e_d = -vrel./vreln;
Fd = 1/2.*b.*vreln.^2.*e_d;

% grid fin control force
gfqs = 1/2.*param.gamma.*pres.*mach.^2.*param.gfS;
Fgfbn = gfqs.*param.dclda(mach).*gf_aoa;
Fgfb = param.gf_point*Fgfbn;
Fgf = C*Fgfb;

% moment arms
rcp = param.rcp_1st;
renge = [-param.d.enge1st; 0; 0];
rc = param.rc1st(m,mox,mfuel);
moic = param.moic1st(m,mox,mfuel);
dmoic = param.dmoic1st(m,mox,mfuel,-Tmax/cstar,0,0);
moicinv = inv(moic);
rcg2cp = rcp - rc;
rcg2eg = renge - rc;

% % % Compute Jacobian for intermediate variables % % %
% dotCmat = [q0 q3 -q2; q1 q2 q3; -q2 q1 -q0; -q3 q0 q1;
%     -q3 q0 q1; q2 -q1 q0; q1 q2 q3; -q0 -q3 q2;
%     q2 -q1 q0; q3 -q0 -q1; q0 q3 -q2; q1 q2 q3]';
dotCmat = [q0 q3 -q2; -q3 q0 q1; q2 -q1 q0;
    q1 q2 q3; q2 -q1 q0; q3 -q0 -q1;
    -q2 q1 -q0; q1 q2 q3; q0 q3 -q2;
    -q3 q0 q1; -q0 -q3 q2; q1 q2 q3];

dotCinvmat = [q0 -q3 q2; q3 q0 -q1; -q2 q1 q0;
    q1 q2 q3; q2 -q1 -q0; q3 q0 -q1; 
    -q2 q1 q0; q1 q2 q3; -q0 q3 -q2;
    -q3 -q0 q1; q0 -q3 q2; q1 q2 q3];

skew = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
RCG2CP = skew(rcg2cp);
RCG2ENG = skew(rcg2eg);

% % Forces % %
% inertial frame - gravity force
Jfgr = m*(3*mu/((r'*r)^2.5)*(r*r') - mu/((r'*r)^1.5)*eye(3));
Jfgv = zeros(3);
Jfgq = zeros(3,4);
Jfgw = zeros(3);
Jfgm = -mu/((r'*r)^1.5)*r;

Jfgy = zeros(3,1);
Jfgp = zeros(3,1);
Jfggfa = zeros(3,4);

% inertial frame - drag force
Jfdr = b/2*(vreln*OMEGA - e_d*(OMEGA*vrel)');
Jfdv = b/2*(-e_d*vrel'-vreln*eye(3));
Jfdq = zeros(3,4);
Jfdw = zeros(3);
Jfdm = zeros(3,1);

Jfdy = zeros(3,1);
Jfdp = zeros(3,1);
Jfdgfa = zeros(3,4);

% body frame - drag force
Jfdbr = C'*Jfdr;
Jfdbv = C'*Jfdv;
Jfdbq = reshape(dotCinvmat*Fd,[3,4]);
Jfdbw = zeros(3);
Jfdbm = zeros(3,1);

Jfdby = zeros(3,1);
Jfdbp = zeros(3,1);
Jfdbgfa = zeros(3,4);

% body frame - thrust
Jftbr = zeros(3);
Jftbv = zeros(3);
Jftbq = zeros(3,4);
Jftbw = zeros(3);
Jftbm = zeros(3,1);

Jftby = Tmax*[-cp*sy; cp*cy; 0];
Jftbp = Tmax*[-cy*sp; -sp*sy; -cp];
Jftbgfa = zeros(3,4);

% inertial frame - thrust
Jftr = zeros(3);
Jftv = zeros(3);
Jftq = reshape(dotCmat*Ftb,[3,4]);
Jftw = zeros(3);
Jftm = zeros(3,1);

Jfty = C*Jftby;
Jftp = C*Jftbp;
Jftgfa = zeros(3,4);

% body frame - grid fin control force
Jfgfbr = zeros(3);
Jfgfbv = zeros(3);
Jfgfbq = zeros(3,4);
Jfgfbw = zeros(3);
Jfgfbm = zeros(3,1);

Jfgfby = zeros(3,1);
Jfgfbp = zeros(3,1);
Jfgfbgfa = param.gf_point*gfqs;

% inertial frame - grid fin control force
Jfgfr = zeros(3);
Jfgfv = zeros(3);
Jfgfq = reshape(dotCmat*Fgfb,[3,4]);
Jfgfw = zeros(3);
Jfgfm = zeros(3,1);

Jfgfy = zeros(3,1);
Jfgfp = zeros(3,1);
Jfgfgfa = C*Jfgfbgfa;

% % Moments % %
% body frame - drag moment
Jmdbr = RCG2CP*Jfdbr;
Jmdbv = RCG2CP*Jfdbv;
Jmdbq = RCG2CP*Jfdbq;
Jmdbw = zeros(3);
Jmdbm = zeros(3,1);

Jmdby = zeros(3,1);
Jmdbp = zeros(3,1);
Jmdbgfa = zeros(3,4);

% body frame - thrust moment
Jmtbr = zeros(3);
Jmtdbv = zeros(3);
Jmtbq = zeros(3,4);
Jmtbw = zeros(3);
Jmtbm = zeros(3,1);

Jmtby = RCG2ENG*Jftby;
Jmtbp = RCG2ENG*Jftbp;
Jmtbgfa = zeros(3,4);

% body frame - grid fin control moment
Jmgfbr = zeros(3);
Jmgfdbv = zeros(3);
Jmgfbq = zeros(3,4);
Jmgfbw = zeros(3);
Jmgfbm = zeros(3,1);

Jmgfby = zeros(3,1);
Jmgfbp = zeros(3,1);
Jmgfbgfa = RCG2ENG*Jfgfbgfa;

% % Angular Momentum Components % %
% I x omega
Jiwr = zeros(3);
Jiwv = zeros(3);
Jiwq = zeros(3,4);
Jiww = moic;
Jiwm = dmoic*omega;

Jiwy = zeros(3,1);
Jiwp = zeros(3,1);
Jiwgfa = zeros(3,4);

% I x omega dot
W = skew(omega);
IW = skew(moic*omega);

Jiwdr = Jmdbr;
Jiwdv = Jmdbv;
Jiwdq = Jmdbq;
Jiwdw = IW-W*Jiww;
Jiwdm = -W*Jiwm;

Jiwdy = Jmtby;
Jiwdp = Jmtbp;
Jiwdgfa = Jmgfbgfa;

% % % Compute state matrix A % % %
% Jacobian - rdot terms
Arr = zeros(3);
Arv = eye(3);
Arq = zeros(3,4);
Arw = zeros(3);
Arm = zeros(3,1);

% Jacobian - vdot terms
Avr = (Jfgr + Jfdr)/m;
Avv = Jfdv/m;
Avq = (Jftq + Jfgfq)/m;
Avw = zeros(3);
Avm = -1/m^2*(Ft+Fd+Fgf);

% Jacobian - qdot terms
Aqr = zeros(4,3);
Aqv = zeros(4,3);
Aqq = [0 -wx -wy -wz; wx 0 wz -wy; wy -wz 0 wx; wz wy -wx 0];
Aqw = [-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0];
Aqm = zeros(4,1);

% Jacobian - wdot terms
Awr = moicinv*Jiwdr;
Awv = moicinv*Jiwdv;
Awq = moicinv*Jiwdq;
Aww = moicinv*Jiwdw;
Awm = moicinv*Jiwdm;

% Jacobian - mdot terms
Amr = zeros(1,3);
Amv = zeros(1,3);
Amq = zeros(1,4);
Amw = zeros(1,3);
Amm = 0;

% % A matrix % %
A = [Arr Arv Arq Arw Arm;
    Avr Avv Avq Avw Avm;
    Aqr Aqv Aqq Aqw Aqm;
    Awr Awv Awq Aww Awm;
    Amr Amv Amq Amw Amm];

% % % Compute input matrix B % % %
% Jacobian - rdot terms
Bry = zeros(3,1);
Brp = zeros(3,1);
Brgfa = zeros(3,4);

% Jacobian - vdot terms
Bvy = Jfty/m;
Bvp = Jftp/m;
Bvgfa = Jfgfgfa/m;

% Jacobian - qdot terms
Bqy = zeros(4,1);
Bqp = zeros(4,1);
Bqgfa = zeros(4,4);

% Jacobian - wdot terms
Bwy = moicinv*Jiwdy;
Bwp = moicinv*Jiwdp;
Bwgfa = moicinv*Jiwdgfa;

% Jacobian - mdot terms
Bmy = 0;
Bmp = 0;
Bmgfa = zeros(1,4);

B = [Bry Brp Brgfa;
    Bvy Bvp Bvgfa;
    Bqy Bqp Bqgfa;
    Bwy Bwp Bwgfa;
    Bmy Bmp Bmgfa];

end