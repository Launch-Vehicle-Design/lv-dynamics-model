%% FUNCTION - Dynamics of a launch vehicle
function [dxdt,interactions] = rocket6dof(t,x,u,param,phase)
% % x represents states as follows nx1
% 1:3 - r, position of the launch vehicle, iner
% 4:6 - v, velocity of the launch vehicle, iner
% 7:10 - q, quaternion from iner frame to body frame
% 11:13 - omega, angular velocity of the body, iner
% 14 - m, mass of the launch vehicle
% 15 - mox, mass of the oxidizer
% 16 - mfuel, mass of the fuel
% 17:18 - phiox and thetaox, body
%   phiox, oxidizer sloshing mass off-b2 pitch anlge
%   thetaox, oxidizer sloshing mass off-b3 roll angle
% 19:20 - phifuel and thetafuel, body
%   phifuel, fuel sloshing mass off-b2 pitch anlge
%   thetafuel, fuel sloshing mass off-b3 roll angle
% 21:22 - dphiox and dthetaox, body
% 23:24 - dphifuel and dthetafuel, body

% % u represents control as follows mx1
% 1:2 - TVC yaw and pitch angles, body
% 3 - "throttle" of main engine
% 4:11 - on-off of the ACS thrusters
% 12:15 - grid fin angles of attacks

% phase indicates operating stages of the vehicle
% 0 - drop coast phase with no main engine firing
% 1 - full stack vehicle with 1st stage burning
% 1.5 - stage separation coast phase with no main engine firing
% 2 - 2nd stage with payload fairing with 2nd stage burning
% 3 - 2nd stage w/o payload fairing with 2nd stage burning

dyn_size = size(x,2);

% decode dynamic states
r = x(1:3,:); v = x(4:6,:); q = x(7:10,:); omega = x(11:13,:);
m = x(14,:); mox = x(15,:); mfuel = x(16,:);
oxang = x(17:18,:); fuelang = x(19:20,:);
doxang = x(21:22,:); dfuelang = x(23:24,:);

% decode controls
tvc_angle = u(1:2,:); throtl = u(3,:);
onoff = u(4:11,:); gf_aoa = u(12:15,:);

% phase indice
ind_phase0 = phase == 0;
ind_phase1 = phase == 1;
ind_phase1p5 = phase == 1.5;
ind_phase2 = phase == 2;
ind_phase3 = phase == 3;
ind_phase01 = ind_phase0 | ind_phase1;
ind_phase23 = ind_phase2 | ind_phase3;

% norm array
rn = vecnorm(r);

% regulate quaternion
q = q./vecnorm(q);
C = EP2C(q);

% atmospheric parameters
h = rn-param.earthR;
atmo_profile = atmo(h);
rho = atmo_profile.rho;
pres = atmo_profile.P;
a = sqrt(param.gamma.*pres./rho);
oxratio = param.OFratio/(1+param.OFratio);
fuelratio = 1-oxratio;

% thrust vectoring pointing
tvc_body = param.tvc_point(tvc_angle(1),tvc_angle(2));
fsfa = vecnorm(tvc_body(2:3))./tvc_body(1);

% aerodynamic drag
atmo_v = param.skewOMEGA*r;
vrel = v - atmo_v;
vreln = vecnorm(vrel);
e_d = -vrel./vreln;
mach = vreln./a;
Cd = 0.5*ones(size(mach));

% components of external forces
Fg = -m.*param.mu./rn.^3.*r;
Fd = 1/2.*rho.*vreln.^2.*Cd*param.S.*e_d;
Facs = C*param.T_acs*param.acs_point*onoff;
mfueldot = vecnorm(Facs)./param.Isp_acs./param.g0;

% % % phase based force computation % % %
% initialization - thrust and grid fin control force
mdot = zeros(size(m)); moxdot = zeros(size(m));
Ftb = zeros(size(Fg)); Ft = zeros(size(Fg));
Fgfbn = zeros(param.gfn,dyn_size); Fgf = zeros(size(Fg)); gfqs = zeros(size(mach));

% % 1st stage % %
% thrust
Ftb(:,ind_phase1) = throtl(:,ind_phase1).*param.maxT_1st.*tvc_body;
mdot(:,ind_phase1) = vecnorm(Ft(:,ind_phase1))./param.Isp1./param.g0 + mfueldot(:,ind_phase1);
% grid fin control force
gfqs(:,ind_phase1) = 1/2.*param.gamma.*pres(:,ind_phase1).*mach(:,ind_phase1).^2.*param.gfS;
Fgfbn(:,ind_phase1) = gfqs(:,ind_phase1).*param.dclda(mach(:,ind_phase1)).*gf_aoa(:,ind_phase1);
Fgf(:,ind_phase1) = C*param.gf_point*Fgfbn(:,ind_phase1);

% % 2nd stage % %
% thrust
Ftb(:,ind_phase23) = throtl(:,ind_phase23).*param.maxT_2nd.*tvc_body;
mpdot = vecnorm(Ft(:,ind_phase23))./param.Isp2./param.g0;
% distribute the ox and fuel mass rate
moxdot(:,ind_phase23) = (oxratio + param.litvc_mdot(fsfa(:,ind_phase23))).*mpdot;
mfueldot(:,ind_phase23) = mfueldot(:,ind_phase23) + fuelratio*mpdot;
mdot(:,ind_phase23) = mfueldot(:,ind_phase23) + moxdot(:,ind_phase23);

% express thrust in iner frame
Ft = C*Ftb;

% Newtonian translational dynamics
atot = (Fg + Fd + Ft + Facs + Fgf)./m;

% Eulerian rotational dynamics
dqdt = zeros([4,dyn_size]);
alphatot = zeros([3,dyn_size]);
rcp = param.rcp_2nd*ones([1,dyn_size]);
renge = [-param.d.enge2nd; 0; 0]*ones([1,dyn_size]);
rcp(:,ind_phase01) = param.rcp_1st*ones([1,sum(ind_phase01)]);
renge(:,ind_phase01) = [-param.d.enge1st; 0; 0]*ones([1,sum(ind_phase01)]);
for i = 1:dyn_size
    % vehicle dimensions in body frame
    if phase(i) == 0
        rc = param.rc1st(param.m0,param.m.ox2nd,param.m.fuel2nd);
        moic = param.moic1st(param.m0,param.m.ox2nd,param.m.fuel2nd);
        dmoic = param.dmoic1st(param.m0,param.m.ox2nd,param.m.fuel2nd,0,0,0);
    elseif phase(i) == 1
        rc = param.rc1st(m(i),mox(i),mfuel(i));
        moic = param.moic1st(m(i),mox(i),mfuel(i));
        dmoic = param.dmoic1st(m(i),mox(i),mfuel(i),mdot(i),moxdot(i),mfueldot(i));
    elseif phase(i) == 1.5
        rc = param.rc2nd_wplf(param.m02,param.m.ox2nd,param.m.fuel2nd);
        moic = param.moic2nd_wplf(param.m02,param.m.ox2nd,param.m.fuel2nd);
        dmoic = param.dmoic2nd_wplf(param.m02,param.m.ox2nd,param.m.fuel2nd,0,0,0);
    elseif phase(i) == 2
        rc = param.rc2nd_wplf(m(i),mox(i),mfuel(i));
        moic = param.moic2nd_wplf(m(i),mox(i),mfuel(i));
        dmoic = param.dmoic2nd_wplf(m(i),mox(i),mfuel(i),mdot(i),moxdot(i),mfueldot(i));
    elseif phase(i) == 3
        rc = param.rc2nd(m(i),mox(i),mfuel(i));
        moic = param.moic2nd(m(i),mox(i),mfuel(i));
        dmoic = param.dmoic2nd(m(i),mox(i),mfuel(i),mdot(i),moxdot(i),mfueldot(i));
    end

    % quaternion rotational kinematics
    w1 = omega(1,i); w2 = omega(2,i); w3 = omega(3,i);
    OMEGA = [0 -w1 -w2 -w3; w1 0 w3 -w2; w2 -w3 0 w1; w3 w2 -w1 0];
    dqdt(:,i) = 1/2*OMEGA*q(:,i);

    % moment arms in body frame
    rcg2cp = rcp(:,i) - rc;
    rcg2eng = renge(:,i) - rc;
    rcg2acs = param.r_acs_mount - rc;
    rcg2gf = param.r_gf_mount - rc;

    % moments in body frame
    Md = cross(rcg2cp, C'*Fd(:,i));
    Mt = cross(rcg2eng, Ftb(:,i));
    Macs = cross(rcg2acs,param.T_acs*param.acs_point)*onoff(:,i);
    Mgf = cross(rcg2gf,Fgfbn(:,i)'.*param.gf_point)*ones(4,1);
    % total external moments expressed at center of mass
    Lc = Md + Mt + Macs + Mgf;

    % euler rotation equation of motion
    alphatot(:,i) = moic\(Lc-dmoic*omega(:,i)-cross(omega(:,i),moic*omega(:,i)));
end

% differential of states
dxdt = zeros(size(x));
dxdt(1:3,:) = v;
dxdt(4:6,:) = atot;
dxdt(7:10,:) = dqdt;
dxdt(11:13,:) = alphatot;
dxdt(14,:) = -mdot;
dxdt(15,:) = -moxdot;
dxdt(16,:) = -mfueldot;
dxdt(17:18,:) = doxang;
dxdt(19:20,:) = dfuelang;
dxdt(21:22,:) = zeros(size(doxang));
dxdt(23:24,:) = zeros(size(dfuelang));

interactions = [Fg; Fd; Ft; Facs; Md; Mt; Macs];
